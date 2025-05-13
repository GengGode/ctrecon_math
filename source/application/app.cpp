#include "debugger.h"
#include <fmt/format.h>
#include <implot.h>
#include <iostream>

// 区分角度和弧度
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>

namespace stdex
{
    struct deg;
    // 弧度类型
    struct rad
    {
        double value;

        // 构造函数
        explicit rad(double val = 0) : value(val) {}

        // 转换为角度
        deg to_deg() const;
    };

    // 角度类型
    struct deg
    {
        double value;

        explicit deg(double val = 0) : value(val) {}

        // 转换为弧度
        rad to_rad() const { return rad(value * std::numbers::pi / 180.0); }
    };
    deg rad::to_deg() const
    {
        return deg(value * 180.0 / std::numbers::pi);
    }

    // 运算符重载（示例：弧度相加）
    rad operator+(const rad& a, const rad& b)
    {
        return rad(a.value + b.value);
    }

    // 运算符重载（示例：角度相加）
    deg operator+(const deg& a, const deg& b)
    {
        return deg(a.value + b.value);
    }

    // 用户定义字面量命名空间
    namespace literals
    {
        // 弧度字面量：1.0_rad
        rad operator"" _rad(long double val)
        {
            return rad(static_cast<double>(val));
        }

        // 角度字面量：90.0_deg
        deg operator"" _deg(long double val)
        {
            return deg(static_cast<double>(val));
        }

        // 支持整数版本（如 5_rad）
        rad operator"" _rad(unsigned long long val)
        {
            return rad(static_cast<double>(val));
        }

        deg operator"" _deg(unsigned long long val)
        {
            return deg(static_cast<double>(val));
        }
    } // namespace literals
} // namespace stdex

// 输出运算符（可选）
std::ostream& operator<<(std::ostream& os, const stdex::rad& r)
{
    os << r.value << " rad";
    return os;
}

std::ostream& operator<<(std::ostream& os, const stdex::deg& d)
{
    os << d.value << " deg";
    return os;
}

using namespace stdex::literals;
//
struct 锥形束正投影
{
    // 1. 该类的构造函数和析构函数
    锥形束正投影() { std::cout << "锥形束正投影对象创建" << std::endl; }
    ~锥形束正投影() { std::cout << "锥形束正投影对象销毁" << std::endl; }

    // 2. 该类的成员函数
    void display() const { std::cout << "锥形束正投影显示" << std::endl; }

    void calculate()
    {
        std::cout << "锥形束正投影计算" << std::endl;
        // 进行计算
        params.计算投影射线束角度宽度内均匀分布的射线线束的相对转动角度(conf);
        for (const auto& angle : params.投影射线束角度宽度内均匀分布的射线线束的相对转动角度)
        {
            std::cout << "投影射线束角度: " << angle << std::endl;
        }
        params.计算投影角度位置(conf);
        // for (const auto& angle : params.所有投影角度位置)
        // {
        //     std::cout << "投影角度: " << angle << std::endl;
        // }
        params.计算投影点位置(conf);
        // for (const auto& point : params.所有投影点)
        // {
        //     std::cout << "投影点: (" << point.real() << ", " << point.imag() << ")" << std::endl;
        // }
        plot_points.clear();
        for (const auto& point : params.所有投影点)
        {
            plot_points.push_back({ point.real(), point.imag() });
        }

        params.计算投影角内射线线束位置(conf);
        fmt::println("投影角内射线线束位置数量: {}", params.所有投影角内射线线束位置.size());
    }

    struct config
    {
        stdex::deg 投影射线角位置的间隔 = 10.0_deg;
        stdex::deg 投影射线角位置的范围 = 60.0_deg;
        stdex::deg 投影射线角位置的起始 = 0.0_deg;
        stdex::deg 投影射线束角度宽度 = 60.0_deg;
        int 投影射线束角度宽度内均匀分布的射线线束的数量 = 6;
        double 锥形束顶点到原点距离 = 100.0;
        double 锥形束的射线线束的长度 = 200.0;
    };
    struct params
    {
        std::vector<stdex::deg> 投影射线束角度宽度内均匀分布的射线线束的相对转动角度;
        void 计算投影射线束角度宽度内均匀分布的射线线束的相对转动角度(config& conf)
        {
            const auto& range = conf.投影射线束角度宽度;
            const auto& count = conf.投影射线束角度宽度内均匀分布的射线线束的数量;
            std::vector<stdex::deg> angls;
            angls.resize(count);
            // 范围等分为count份，然后取其中点
            for (int i = 0; i < count; ++i)
            {
                double angle = range.value / count * (i + 0.5) - range.value / 2;
                angls[i] = stdex::deg(angle);
            }
            投影射线束角度宽度内均匀分布的射线线束的相对转动角度 = angls;
        }
        int 投影角度数量 = 0;
        void 计算投影角度数量(config& conf) { 投影角度数量 = static_cast<int>(conf.投影射线角位置的范围.value / conf.投影射线角位置的间隔.value); }
        std::vector<stdex::deg> 所有投影角度位置;
        void 计算投影角度位置(config& conf)
        {
            所有投影角度位置.clear();
            for (double angle = conf.投影射线角位置的起始.value; angle < conf.投影射线角位置的范围.value; angle += conf.投影射线角位置的间隔.value)
            {
                所有投影角度位置.push_back(stdex::deg(angle));
            }
            投影角度数量 = 所有投影角度位置.size();
        }
        std::vector<std::complex<double>> 所有投影点;
        void 计算投影点位置(config& conf)
        {
            所有投影点.clear();
            for (const auto& angle : 所有投影角度位置)
            {
                所有投影点.push_back(std::polar(conf.锥形束顶点到原点距离, angle.to_rad().value));
            }
        }

        // using ImPlotLine = ImPlotRect;
        std::vector<std::pair<ImPlotPoint, ImPlotPoint>> 所有投影角内射线线束位置;
        void 计算投影角内射线线束位置(config& conf)
        {
            // 每个投影点作为起点，分别各分布的线束的末端为止点
            所有投影角内射线线束位置.clear();
            for (const auto& point : 所有投影点)
            {
                for (const auto& angle : 投影射线束角度宽度内均匀分布的射线线束的相对转动角度)
                {
                    auto arg = std::arg(-point);
                    // 计算每个投影点的线束位置
                    double x = point.real() + conf.锥形束的射线线束的长度 * std::cos(arg + angle.to_rad().value);
                    double y = point.imag() + conf.锥形束的射线线束的长度 * std::sin(arg + angle.to_rad().value);
                    所有投影角内射线线束位置.push_back({ { point.real(), point.imag() }, { x, y } });
                    fmt::println("投影点: ({}, {}), 线束位置: ({}, {})", point.real(), point.imag(), x, y);
                }
            }
        }
    };

    config conf;
    params params;
    std::vector<ImPlotPoint> plot_points;
};

struct plot_render : irender
{
    plot_render(锥形束正投影& projection) : projection(projection) {}
    锥形束正投影& projection;
    void render() override
    {
        ImPlot::BeginPlot("Plot Var Pool", { -1, -1 });
        // for (const auto& angle : projection.params.所有投影角度位置)
        // {
        //     ImPlot::PlotLine("投影角度", &angle.value, 1);
        // }
        // for (const auto& point : projection.plot_points)
        //{
        //    ImPlot::PlotLine("投影点", &point.real(), 1);
        //}
        // 绘制投影点
        ImPlot::PlotScatter("投影点", &projection.plot_points[0].x, &projection.plot_points[0].y, projection.plot_points.size(), 0, 0, sizeof(ImPlotPoint));

        // 绘制投影角内射线线束位置
        for (const auto& line : projection.params.所有投影角内射线线束位置)
        {
            double xs[2] = { line.first.x, line.second.x };
            double ys[2] = { line.first.y, line.second.y };
            ImPlot::PlotLine("投影角内射线线束位置", xs, ys, 2);
        }

        ImPlot::EndPlot();
    }
};

int main()
{
    system("chcp 65001");

    锥形束正投影 projection;
    projection.display();
    projection.calculate();

    debugger d(std::make_shared<plot_render>(projection));
    return d.execute();
}
