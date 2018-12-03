#include "Scanner.h"

using namespace scanner;

static bool sortByArea(const std::vector<cv::Point> &v1, const std::vector<cv::Point> &v2) {
    double v1Area = fabs(contourArea(cv::Mat(v1)));
    double v2Area = fabs(contourArea(cv::Mat(v2)));
    return v1Area > v2Area;
}

//计算四个点围成四边形的面积，排序比较
static bool sortArea( std::vector<cv::Point> &p1, std::vector<cv::Point> &p2){
    cv::Point leftTop = p1[0];
    cv::Point rightTop = p1[2];
    cv::Point rightBottom = p1[1];
    cv::Point leftBottom = p1[3];
    double s1 = fabs( 1/2*(leftTop.x*rightTop.y+rightBottom.x*rightBottom.y+leftBottom.x*leftTop.y-leftTop.x*rightBottom.y-rightTop.x*leftBottom.y-leftBottom.x*leftTop.y));
    double s2 = 1/2*(leftTop.x*rightBottom.y+rightBottom.x*leftBottom.y+leftBottom.x*leftTop.y-leftTop.x*leftBottom.y-rightBottom.x*leftTop.y-leftBottom.x*rightBottom.y);
    double sp1 = s1+s2;

    cv::Point leftTop2 = p2[0];
    cv::Point rightTop2 = p2[1];
    cv::Point rightBottom2 = p2[2];
    cv::Point leftBottom2 = p2[3];
    double s3 = 1/2*(leftTop2.x*rightTop2.y+rightBottom2.x*rightBottom2.y+leftBottom2.x*leftTop2.y-leftTop2.x*rightBottom2.y-rightTop2.x*leftBottom2.y-leftBottom2.x*leftTop2.y);
    double s4 = 1/2*(leftTop2.x*rightBottom2.y+rightBottom2.x*leftBottom2.y+leftBottom2.x*leftTop2.y-leftTop2.x*leftBottom2.y-rightBottom2.x*leftTop2.y-leftBottom2.x*rightBottom2.y);
    double sp2 = s3+s4;

    return sp1 > sp2;
}



Scanner::Scanner(cv::Mat &bitmap) {
    srcBitmap = bitmap;
}

Scanner::~Scanner() {
}

std::vector<cv::Point> Scanner::scanPoint() {
    //缩小图片尺寸
    cv::Mat image = resizeImage();
    //预处理图片
    cv::Mat scanImage = preprocessImage(image);
    std::vector<std::vector<cv::Point>> contours;
    //提取边框
    findContours(scanImage, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
    std::vector<cv::Point2f> conerpoint;
    cv::Mat empty=cv::Mat();
    //empty.setTo(255);
    cv::goodFeaturesToTrack(scanImage,conerpoint,200,0.01,25,empty,2,false,0.04);//mask需要设置255
    //按面积排序
    sort(contours.begin(), contours.end(), sortByArea);
    std::vector<cv::Point> result;
    if (contours.size() > 0) {
        std::vector<cv::Point> contour = contours[0];//youwenti
        double arc = arcLength(contour, true);

        //找距离最远的两个点
        cv::Point lp1;
        cv::Point lp2;
        bool first = true;
        double temp;
        for(int i = 0;i<contour.size()-1;i++){
            for (int j = i+1; j<contour.size();j++){
                double len = getdis(contour[i],contour[j]);
                if(len>temp){
                    temp = len;
                    lp1.x=contour[i].x;
                    lp1.y=contour[i].y;
                    lp2.x=contour[j].x;
                    lp2.y=contour[j].y;
                }

            }

        }
        //筛选轮廓点重复部分
        //checkpoint(contour,contour,lass);
        //获取角点
        //std::vector<cv::Point> outDP;
        //获取角点（有效）
        //shaPoint(contour,outDP);
        std::vector<cv::Point> finalp;
        filpoint(conerpoint,contour,finalp);

        //将角点进行排列组合
        int num = 2;
        int result1[num];
        std::vector<cv::Point> combinpoints;
        combinepoint(finalp,0,result1,num,num,finalp.size(),combinpoints);

        int n = combinpoints.size()/2;
        int m = 2;
        std::vector<std::vector<cv::Point> > compoints(n,std::vector<cv::Point>(0));
        std::vector<std::vector<cv::Point> > parepoints(n,std::vector<cv::Point>(0));
        for(int i=0;i<n;i++){
            compoints[i].push_back(lp1);
            compoints[i].push_back(lp2);
            for(int j = 0;j<m;j++){

                compoints[i].push_back(combinpoints[i*2+j]);
            }
        }
        int all=0;

        for(int i = 0;i<compoints.size();i++){
        if(compareDis(compoints,i) == true){
                for(int j = 0;j<4;j++){
                    parepoints[all].push_back(compoints[i][j]);
                }
            all++;
        }
        }
       cv::Point fp1,fp2,fp3,fp4;
        double tmp = 1.11;
        for(int i = 0; i<all; i++){
            std::vector<cv::Point> &p1 = parepoints[i];
            double s1 = fabs(1/2*(p1[0].x*p1[2].y+p1[1].x*p1[1].y+p1[3].x*p1[0].y-p1[0].x*p1[1].y-p1[2].x*p1[3].y-p1[3].x*p1[0].y));
            double s2 = fabs(1/2*(p1[0].x*p1[1].y+p1[1].x*p1[3].y+p1[3].x*p1[0].y-p1[0].x*p1[3].y-p1[1].x*p1[2].y-p1[3].x*p1[1].y));
            double sp1 = s1+s2;
           //double area1 = getArea(parepoints[i]);
            if(sp1>tmp){
                tmp = sp1;
                fp1.x=parepoints[i][0].x;
                fp2.x=parepoints[i][1].x;
                fp3.x=parepoints[i][2].x;
                fp4.x=parepoints[i][3].x;
                fp1.y=parepoints[i][0].y;
                fp2.y=parepoints[i][1].y;
                fp3.y=parepoints[i][2].y;
                fp4.y=parepoints[i][3].y;

            }
        }



        sort(parepoints.begin(),parepoints.end(),sortArea);

        //approxPolyDP(cv::Mat(contour), outDP, 0.02 * arc, true);
        //compareDis(conerpoint,ss);


            //筛选去除相近的点
            std::vector<cv::Point> selectedPoints = selectPoints(finalp, 1);
            if (selectedPoints.size() != 4) {
                //如果筛选出来之后不是四边形，那么使用最小矩形包裹
                cv::RotatedRect rect = minAreaRect(contour);
                //cv::RotatedRect rect = minAreaRect(finalp);
                cv::Point2f p[4];
                rect.points(p);
                result.push_back(p[0]);
                result.push_back(p[1]);
                result.push_back(p[2]);
                result.push_back(p[3]);
            } else {
                result = selectedPoints;
            }
            for (cv::Point &p : result) {
                p.x *= resizeScale;
                p.y *= resizeScale;
        }
    }
    // 按左上，右上，右下，左下排序，方便下面使用凸包进行计算
    std::vector<cv::Point> r = sortPointClockwise(result);
    //如果不是凸包，则用rect进行包裹，保证输出的是凸包
//    if(gimp_transform_polygon_is_convex(r)<=0){
//        cv::RotatedRect rect = minAreaRect(result);
//        cv::Point2f p[4];
//        rect.points(p);
//        result.clear();
//        result.push_back(p[0]);
//        result.push_back(p[1]);
//        result.push_back(p[2]);
//        result.push_back(p[3]);
//    }
//    return r;
}



cv::Mat Scanner::resizeImage() {
    int width = srcBitmap.cols;
    int height = srcBitmap.rows;
    int maxSize = width > height ? width : height;
    if (maxSize > resizeThreshold) {
        resizeScale = 1.0f * maxSize / resizeThreshold;
        width = static_cast<int>(width / resizeScale);
        height = static_cast<int>(height / resizeScale);
        cv::Size size(width, height);
        cv::Mat resizedBitmap(size, CV_8UC3);
        resize(srcBitmap, resizedBitmap, size);
        return resizedBitmap;
    }
    return srcBitmap;
}



cv::Mat Scanner::preprocessImage(cv::Mat &image) {
    cv::Mat grayMat;
    cvtColor(image, grayMat, CV_BGR2GRAY);
    cv::Mat blurMat;
    GaussianBlur(grayMat, blurMat, cv::Size(5, 5), 0);
    cv::Mat cannyMat;
    Canny(blurMat, cannyMat, 10, 20 );
    cv::Mat thresholdMat;
    threshold(cannyMat, thresholdMat, 0, 255, CV_THRESH_OTSU);
    return thresholdMat;
}

std::vector<cv::Point> Scanner::selectPoints(std::vector<cv::Point> points, int selectTimes) {
    if (points.size() > 4) {
        double arc = arcLength(points, true);
        std::vector<cv::Point>::iterator itor = points.begin();
        while (itor != points.end()) {
            if (points.size() == 4) {
                return points;
            }
            cv::Point &p = *itor;
            if (itor != points.begin()) {
                cv::Point &lastP = *(itor - 1);
                double pointLength = sqrt(pow((p.x - lastP.x), 2) + pow((p.y - lastP.y), 2));
                if (pointLength < arc * 0.1 * selectTimes && points.size() > 4) {
                    itor = points.erase(itor);
                    continue;
                }
            }
            itor++;
        }
        if (points.size() > 4) {
            return selectPoints(points, selectTimes + 1);
        }
    }
    return points;
}

//修改了部分内容，之前的点太多了
//bool Scanner::checkpoint(std::vector<cv::Point> po1, std::vector<cv::Point> po2,std::vector<cv::Point> &getPoints){
//
//
//    for (int i = 0;i<po1.size();i++){
//        for(int j = i+1;j<po1.size();j++){
//            if(po1[i].x-po1[j].x==0 && po1[i].y-po1[j].y){
//                for(int k = j+1;k<po1.size();k++){
//                    po1[k-1]==po1[k];
//                }
//                po1.size()-1;
//                j--;
//            }
//        }
//
//        }
//
//    return true;
//}



bool Scanner::filpoint(std::vector<cv::Point2f> points1, std::vector<cv::Point> points2,std::vector<cv::Point> &getPoint) {
    //int i = 0, j = 0,m=0;

    for(int i=0;i<points1.size();i++){
        for (int j=0;j<points2.size();j++){
            if (fabs(points1[i].x - points2[j].x) <1.5 && fabs(points1[i].y - points2[j].y < 1.5)) {
                getPoint.push_back(points2[j]);//修改返回的数据为int型
            }
            else{
                continue;
            }
        }
    }
//    while (i < points1.size() && j < points2.size()) {
//        if (fabs(points1[i].x - points2[j].x) <2 && fabs(points1[i].y - points2[j].y < 2)) {
//            /getPoint.push_back(points2[j]);//修改返回的数据为int型
//            i++;
//           j = 0;
//        }
//        else {
//            j++;
//        }
//    }
    return true;
}


std::vector<cv::Point> Scanner::sortPointClockwise(std::vector<cv::Point> points) {
    if (points.size() != 4) {
        return points;
    }

    cv::Point unFoundPoint;
    std::vector<cv::Point> result = {unFoundPoint, unFoundPoint, unFoundPoint, unFoundPoint};

    long minDistance = -1;
    for (cv::Point &point : points) {
        long distance = point.x * point.x + point.y * point.y;
        if (minDistance == -1 || distance < minDistance) {
            result[0] = point;
            minDistance = distance;
        }
    }
    if (result[0] != unFoundPoint) {
        cv::Point &leftTop = result[0];
        points.erase(std::remove(points.begin(), points.end(), leftTop));
        if ((pointSideLine(leftTop, points[0], points[1]) *
             pointSideLine(leftTop, points[0], points[2])) < 0) {
            result[2] = points[0];
        } else if ((pointSideLine(leftTop, points[1], points[0]) *
                    pointSideLine(leftTop, points[1], points[2])) < 0) {
            result[2] = points[1];
        } else if ((pointSideLine(leftTop, points[2], points[0]) *
                    pointSideLine(leftTop, points[2], points[1])) < 0) {
            result[2] = points[2];
        }
    }
    if (result[0] != unFoundPoint && result[2] != unFoundPoint) {
        cv::Point &leftTop = result[0];
        cv::Point &rightBottom = result[2];
        points.erase(std::remove(points.begin(), points.end(), rightBottom));
        if (pointSideLine(leftTop, rightBottom, points[0]) > 0) {
            result[1] = points[0];
            result[3] = points[1];
        } else {
            result[1] = points[1];
            result[3] = points[0];
        }
    }

    if (result[0] != unFoundPoint && result[1] != unFoundPoint && result[2] != unFoundPoint &&
        result[3] != unFoundPoint) {
        return result;
    }

    return points;
}

long long Scanner::pointSideLine(cv::Point &lineP1, cv::Point &lineP2, cv::Point &point) {
    long x1 = lineP1.x;
    long y1 = lineP1.y;
    long x2 = lineP2.x;
    long y2 = lineP2.y;
    long x = point.x;
    long y = point.y;
    return (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1);
}

//判断四边形是否为凸包
//int Scanner::gimp_transform_polygon_is_convex(std::vector<cv::Point> vector) {
//    double z1, z2, z3, z4;
//    vector = sortPointClockwise(vector);
//    cv::Point &leftTop = vector[0];
//    cv::Point &rightTop = vector[1];
//    cv::Point &rightBottom = vector[2];
//    cv::Point &leftBottom = vector[3];
//
//
//    z1 = ((rightTop.x - leftTop.x) * (rightBottom.y - leftTop.y) -
//          (rightBottom.x - leftTop.x) * (rightTop.y - leftTop.y));
//    z2 = ((rightBottom.x - leftTop.x) * (leftBottom.y - leftTop.y) -
//          (leftBottom.x - leftTop.x) * (rightBottom.y - leftTop.y));
//    z3 = ((rightBottom.x - rightTop.x) * (leftBottom.y - rightTop.y) -
//          (leftBottom.x - rightTop.x) * (rightBottom.y - rightTop.y));
//    z4 = ((leftBottom.x - rightTop.x) * (leftTop.y - rightTop.y) -
//          (leftTop.x - rightTop.x) * (leftBottom.y - rightTop.y));
//
//    return (z1 * z2 > 0) && (z3 * z4 > 0);
//}




//角点排列组合
bool Scanner::combinepoint(std::vector<cv::Point> points,int start,int* result,int count,const int num,const int length,std::vector<cv::Point> &combine){
    int i =0;
    for(i=start;i<length+1-count;i++){
        result[count-1]=i;
        if(count-1==0){
            int j;
            int n;
            for(j=num-1;j>=0;j--){
                combine.push_back(points[result[j]]);//函数传入出现问题。如何传入此类型函数？？？？

            }


        }
        else{
            combinepoint(points,i+1,result,count-1, num, length, combine);
        }
    }
    return true;
}

double Scanner::calArea(std::vector<cv::Point> vector){
    cv::Point &leftTop = vector[0];
    cv::Point &rightTop = vector[1];
    cv::Point &rightBottom = vector[2];
    cv::Point &leftBottom = vector[3];
    double s1 = 1/2*(leftTop.x*rightTop.y+rightBottom.x*rightBottom.y+leftBottom.x*leftTop.y-leftTop.x*rightBottom.y-rightTop.x*leftBottom.y-leftBottom.x*leftTop.y);
    double s2 = 1/2*(leftTop.x*rightBottom.y+rightBottom.x*leftBottom.y+leftBottom.x*leftTop.y-leftTop.x*leftBottom.y-rightBottom.x*leftTop.y-leftBottom.x*rightBottom.y);

    return (s1+s2);
}




double Scanner::getdis( cv::Point p1,  cv::Point p2){

    return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.x-p2.y)*(p1.y-p2.y));
}

double Scanner::getArea(std::vector<cv::Point> &p1){
    cv::Point leftTop = p1[2];
    cv::Point rightTop = p1[2];
    cv::Point rightBottom = p1[1];
    cv::Point leftBottom = p1[3];
    double s1 = fabs(1/2*(p1[0].x*p1[2].y+p1[1].x*p1[1].y+p1[3].x*p1[0].y-p1[0].x*p1[1].y-p1[2].x*p1[3].y-p1[3].x*p1[0].y));
    double s2 = fabs(1/2*(p1[0].x*p1[1].y+p1[1].x*p1[3].y+p1[3].x*p1[0].y-p1[0].x*p1[3].y-p1[1].x*p1[2].y-p1[3].x*p1[1].y));
    double sp1 = s1+s2;

    return sp1;

}



bool Scanner::compareDis(std::vector<std::vector<cv::Point>> points,int num) {
    int i = 0, j = i + 1, m = j + 1, n = m + 1;
    double b1=0.45;
    double b2=2.222;
    double b3=0.667;
    double b4=1.5;


    double a = fabs(sqrt(pow((points[num][i].x - points[num][n].x), 2) +
                                    pow((points[num][i].y - points[num][n].y), 2)));
    double b = fabs(sqrt(pow((points[num][i].x - points[num][m].x), 2) +
                                    pow((points[num][i].y - points[num][m].y), 2)));
    double c = fabs(sqrt(pow((points[num][j].x - points[num][m].x), 2) +
                                    pow((points[num][j].y - points[num][m].y), 2)));
    double d = fabs(sqrt(pow((points[num][j].x - points[num][n].x), 2) +
                                    pow((points[num][j].y - points[num][n].y), 2)));
    if (b1 < a / b && a/b<b2 && b1 < b / c && b / c< b2 && b1< c / d &&c / d< b2 && b1 < a / d &&a / d< b2 && b3 < a / c &&a / c< b4 && b3 < b / d&&b / d < b4 ) {

        return true;
    }else{
        return false;
    }


}

