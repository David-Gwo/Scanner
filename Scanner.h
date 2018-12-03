#ifndef MONITOR_SCANNER_H
#define MONITOR_SCANNER_H


#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>

namespace scanner{

    class Scanner {
    public:
        int resizeThreshold = 500;

        Scanner(cv::Mat& bitmap);
        virtual ~Scanner();
        std::vector<cv::Point> scanPoint();
    private:
        cv::Mat srcBitmap;
        float resizeScale = 1.0f;

        cv::Mat resizeImage();

        cv::Mat preprocessImage(cv::Mat& image);

        std::vector<cv::Point2f> getConerPoint(cv::Mat &image);

        std::vector<cv::Point> selectPoints(std::vector<cv::Point> points, int selectTimes);


        std::vector<cv::Point> sortPointClockwise(std::vector<cv::Point> vector);

        bool shaPoint(std::vector<cv::Point> vector, std::vector<cv::Point> &result);

        bool checkpoint(std::vector<cv::Point> po1,std::vector<cv::Point> po2, std::vector<cv::Point> &getPoints);

        bool compareDis(std::vector<cv::Point2f> points, std::vector<cv::Point> &repoint);

        bool filpoint(std::vector<cv::Point2f> points1, std::vector<cv::Point> points2,std::vector<cv::Point> &getPoint);

        int gimp_transform_polygon_is_convex(std::vector<cv::Point> vector);

        long long pointSideLine(cv::Point& lineP1, cv::Point& lineP2, cv::Point& point);
    };



}


#endif //MONITOR_SCANNER_H
