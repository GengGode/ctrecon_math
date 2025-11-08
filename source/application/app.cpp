#include "debugger.h"
#include <fmt/format.h>
#include <implot.h>
#include <implot3d.h>
#include <iostream>
#include <vector>
#include <cstdio>

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

// 三维点结构
struct Point3D
{
    double x, y, z;
    Point3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
};

// 螺旋锥束CT扫描器（基于论文）
struct 螺旋锥束扫描器
{
    // 论文中的几何参数
    struct 几何参数
    {
        double RF = 570.0;  // 焦点到旋转中心的距离 (mm)
        double RD = 435.0;  // 探测器到旋转中心的距离 (mm)
        double d = 10.0;    // 每360°旋转的床进增量 (mm)
        double RM = 250.0;  // 测量场半径 (mm)
    };

    几何参数 参数;

    // 计算螺旋焦点轨迹 s(a)
    // 公式(1): s(a) = [RF*sin(a), -RF*cos(a), d*a/(2π)]
    Point3D 计算焦点位置(double a) const
    {
        const double two_pi = 2.0 * std::numbers::pi;
        return Point3D(
            参数.RF * std::sin(a),
            -参数.RF * std::cos(a),
            参数.d * a / two_pi
        );
    }

    // 计算探测器平面上的点 r(u,v)
    // 公式(2a): D: r(u,v) = [-RD*sin(a), RD*cos(a), d*a/(2π)] + u*[cos(a), sin(a), 0] + v*[0, 0, 1]
    Point3D 计算探测器点(double a, double u, double v) const
    {
        const double two_pi = 2.0 * std::numbers::pi;
        const double sin_a = std::sin(a);
        const double cos_a = std::cos(a);
        
        return Point3D(
            -参数.RD * sin_a + u * cos_a,
            参数.RD * cos_a + u * sin_a,
            参数.d * a / two_pi + v
        );
    }

    // 计算从焦点到探测器的射线
    struct 射线
    {
        Point3D 起点;  // 焦点位置
        Point3D 终点;  // 探测器点
    };

    // 生成指定角度位置的射线束
    std::vector<射线> 生成射线束(double a, double u_min, double u_max, int u_count, double v_min, double v_max, int v_count) const
    {
        std::vector<射线> rays;
        rays.reserve(u_count * v_count);

        Point3D focus = 计算焦点位置(a);

        for (int i = 0; i < u_count; ++i)
        {
            double u = u_min + (u_max - u_min) * i / (u_count - 1);
            for (int j = 0; j < v_count; ++j)
            {
                double v = v_min + (v_max - v_min) * j / (v_count - 1);
                Point3D detector = 计算探测器点(a, u, v);
                rays.push_back({ focus, detector });
            }
        }

        return rays;
    }

    // 生成螺旋轨迹点
    std::vector<Point3D> 生成螺旋轨迹(double a_start, double a_end, int count) const
    {
        std::vector<Point3D> trajectory;
        trajectory.reserve(count);

        for (int i = 0; i < count; ++i)
        {
            double a = a_start + (a_end - a_start) * i / (count - 1);
            trajectory.push_back(计算焦点位置(a));
        }

        return trajectory;
    }

    // 计算螺旋线在角度a处的切线方向（单位向量）
    // s'(a) = [RF*cos(a), RF*sin(a), d/(2π)]
    Point3D 计算切线方向(double a) const
    {
        const double two_pi = 2.0 * std::numbers::pi;
        const double cos_a = std::cos(a);
        const double sin_a = std::sin(a);
        
        Point3D tangent(
            参数.RF * cos_a,
            参数.RF * sin_a,
            参数.d / two_pi
        );
        
        // 归一化
        double len = std::sqrt(tangent.x * tangent.x + tangent.y * tangent.y + tangent.z * tangent.z);
        if (len > 1e-10)
        {
            tangent.x /= len;
            tangent.y /= len;
            tangent.z /= len;
        }
        
        return tangent;
    }

    // 计算旋转中心位置（随射线源z坐标移动，保持水平齐平）
    Point3D 计算旋转中心(double a) const
    {
        Point3D focus = 计算焦点位置(a);
        // 中心在z轴上，与射线源保持同一水平高度
        return Point3D(0.0, 0.0, focus.z);
    }

    // 计算从焦点指向中心的方向（径向方向，单位向量）
    Point3D 计算径向方向(double a) const
    {
        Point3D focus = 计算焦点位置(a);
        Point3D center = 计算旋转中心(a);
        
        // 从焦点指向旋转中心
        Point3D radial(center.x - focus.x, center.y - focus.y, center.z - focus.z);
        
        // 归一化
        double len = std::sqrt(radial.x * radial.x + radial.y * radial.y + radial.z * radial.z);
        if (len > 1e-10)
        {
            radial.x /= len;
            radial.y /= len;
            radial.z /= len;
        }
        
        return radial;
    }

    // 生成切线平面上的点
    // 平面由两个方向定义：切线方向和径向方向
    // 平面大小由size参数控制
    std::vector<Point3D> 生成切线平面(double a, double size = 100.0, int resolution = 20) const
    {
        Point3D focus = 计算焦点位置(a);
        Point3D tangent = 计算切线方向(a);
        Point3D radial = 计算径向方向(a);
        
        std::vector<Point3D> plane_points;
        plane_points.reserve(resolution * resolution);
        
        // 生成平面网格点
        for (int i = 0; i < resolution; ++i)
        {
            for (int j = 0; j < resolution; ++j)
            {
                // 在[-size, size]范围内生成点
                double u = -size + (2.0 * size) * i / (resolution - 1);
                double v = -size + (2.0 * size) * j / (resolution - 1);
                
                // 平面上的点 = 焦点 + u*切线方向 + v*径向方向
                Point3D pt(
                    focus.x + u * tangent.x + v * radial.x,
                    focus.y + u * tangent.y + v * radial.y,
                    focus.z + u * tangent.z + v * radial.z
                );
                
                plane_points.push_back(pt);
            }
        }
        
        return plane_points;
    }
};

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

// 3D螺旋轨迹可视化渲染器
struct spiral_3d_render : irender
{
    spiral_3d_render(螺旋锥束扫描器& scanner) : scanner(scanner) 
    {
        // 初始化参数
        当前角度 = 0.0;
        角度步长 = 0.1;
        u范围 = 200.0;  // 探测器u方向范围 (mm)
        v范围 = 50.0;   // 探测器v方向范围 (mm)
        u采样数 = 10;
        v采样数 = 5;
        轨迹圈数 = 2.0;
        动画开启 = false;
        动画速度 = 1.0;  // 角度/秒
        上次更新时间 = 0.0;
        
        // 生成初始数据
        更新数据();
    }

    螺旋锥束扫描器& scanner;
    
    // 控制参数
    double 当前角度;
    double 角度步长;
    double u范围;
    double v范围;
    int u采样数;
    int v采样数;
    double 轨迹圈数;
    
    // 动画参数
    bool 动画开启;
    double 动画速度;  // 弧度/秒
    double 上次更新时间;
    
    // 数据缓存
    std::vector<Point3D> 螺旋轨迹点;
    std::vector<螺旋锥束扫描器::射线> 当前射线束;
    std::vector<Point3D> 切线平面点;
    bool 显示切线平面 = true;
    double 切线平面大小 = 100.0;
    
    void 更新数据()
    {
        // 生成螺旋轨迹（多圈）
        const double two_pi = 2.0 * std::numbers::pi;
        double a_start = 0.0;
        double a_end = 轨迹圈数 * two_pi;
        int 轨迹点数 = static_cast<int>(轨迹圈数 * 360);
        螺旋轨迹点 = scanner.生成螺旋轨迹(a_start, a_end, 轨迹点数);
        
        // 生成当前角度的射线束
        当前射线束 = scanner.生成射线束(
            当前角度,
            -u范围, u范围, u采样数,
            -v范围, v范围, v采样数
        );
        
        // 生成切线平面
        if (显示切线平面)
        {
            切线平面点 = scanner.生成切线平面(当前角度, 切线平面大小, 20);
        }
    }

    void render() override
    {
        ImGui::Begin("螺旋锥束CT扫描器可视化");
        
        // 自动动画更新
        if (动画开启)
        {
            double 当前时间 = ImGui::GetTime();
            if (上次更新时间 > 0.0)
            {
                double 时间差 = 当前时间 - 上次更新时间;
                当前角度 += 动画速度 * 时间差;
                
                const double two_pi = 2.0 * std::numbers::pi;
                double 最大角度 = 轨迹圈数 * two_pi;
                if (当前角度 >= 最大角度)
                {
                    当前角度 = std::fmod(当前角度, 最大角度);
                }
                
                更新数据();
            }
            上次更新时间 = 当前时间;
        }
        else
        {
            上次更新时间 = 0.0;
        }
        
        // 控制面板
        if (ImGui::CollapsingHeader("控制参数"))
        {
            bool 需要更新 = false;
            
            // 动画控制
            if (ImGui::Checkbox("自动动画", &动画开启))
            {
                if (动画开启)
                {
                    上次更新时间 = ImGui::GetTime();
                }
            }
            
            if (动画开启)
            {
                float 动画速度_f = static_cast<float>(动画速度);
                if (ImGui::SliderFloat("动画速度 (rad/s)", &动画速度_f, 0.1f, 10.0f))
                {
                    动画速度 = 动画速度_f;
                }
            }
            
            ImGui::Separator();
            
            float 当前角度_f = static_cast<float>(当前角度);
            if (ImGui::SliderFloat("当前角度 (rad)", &当前角度_f, 0.0f, static_cast<float>(4.0 * std::numbers::pi)))
            {
                当前角度 = 当前角度_f;
                需要更新 = true;
            }
            
            float 角度步长_f = static_cast<float>(角度步长);
            if (ImGui::SliderFloat("角度步长", &角度步长_f, 0.01f, 1.0f))
            {
                角度步长 = 角度步长_f;
                需要更新 = true;
            }
            
            float u范围_f = static_cast<float>(u范围);
            if (ImGui::SliderFloat("探测器U范围 (mm)", &u范围_f, 50.0f, 500.0f))
            {
                u范围 = u范围_f;
                需要更新 = true;
            }
            
            float v范围_f = static_cast<float>(v范围);
            if (ImGui::SliderFloat("探测器V范围 (mm)", &v范围_f, 10.0f, 200.0f))
            {
                v范围 = v范围_f;
                需要更新 = true;
            }
            
            if (ImGui::SliderInt("U采样数", &u采样数, 3, 20))
                需要更新 = true;
            
            if (ImGui::SliderInt("V采样数", &v采样数, 2, 10))
                需要更新 = true;
            
            float 轨迹圈数_f = static_cast<float>(轨迹圈数);
            if (ImGui::SliderFloat("轨迹圈数", &轨迹圈数_f, 0.5f, 5.0f))
            {
                轨迹圈数 = 轨迹圈数_f;
                需要更新 = true;
            }
            
            if (ImGui::Button("更新数据") || 需要更新)
                更新数据();
            
            if (ImGui::Button("重置角度"))
            {
                当前角度 = 0.0;
                更新数据();
            }
            
            ImGui::Separator();
            
            // 切线平面控制
            if (ImGui::Checkbox("显示切线平面", &显示切线平面))
            {
                更新数据();
            }
            
            if (显示切线平面)
            {
                float 平面大小_f = static_cast<float>(切线平面大小);
                if (ImGui::SliderFloat("切线平面大小 (mm)", &平面大小_f, 50.0f, 300.0f))
                {
                    切线平面大小 = 平面大小_f;
                    更新数据();
                }
            }
        }
        
        // 几何参数显示和调整
        if (ImGui::CollapsingHeader("几何参数"))
        {
            bool 参数已更改 = false;
            
            // 螺距参数（pitch = d / S，其中S是层厚）
            static double 层厚 = 1.0;  // 层厚 (mm)，默认1mm
            static double 螺距 = scanner.参数.d / 层厚;  // 初始螺距
            
            float 层厚_f = static_cast<float>(层厚);
            if (ImGui::SliderFloat("层厚 S (mm)", &层厚_f, 0.1f, 10.0f))
            {
                层厚 = 层厚_f;
                // 根据螺距和层厚更新床进增量
                scanner.参数.d = 螺距 * 层厚;
                参数已更改 = true;
            }
            
            float 螺距_f = static_cast<float>(螺距);
            if (ImGui::SliderFloat("螺距 (pitch)", &螺距_f, 0.1f, 5.0f))
            {
                螺距 = 螺距_f;
                // 根据螺距和层厚更新床进增量
                scanner.参数.d = 螺距 * 层厚;
                参数已更改 = true;
            }
            
            // 直接调整床进增量（会自动更新螺距）
            float d_f = static_cast<float>(scanner.参数.d);
            if (ImGui::SliderFloat("床进增量 d (mm/360°)", &d_f, 0.1f, 50.0f))
            {
                scanner.参数.d = d_f;
                // 根据床进增量和层厚更新螺距
                if (层厚 > 1e-6)
                {
                    螺距 = scanner.参数.d / 层厚;
                }
                参数已更改 = true;
            }
            
            ImGui::Separator();
            
            // 显示其他固定参数
            ImGui::Text("RF (焦点距离): %.1f mm", scanner.参数.RF);
            ImGui::Text("RD (探测器距离): %.1f mm", scanner.参数.RD);
            ImGui::Text("RM (测量场半径): %.1f mm", scanner.参数.RM);
            
            // 显示当前计算的螺距值
            ImGui::Text("当前螺距: %.3f (d=%.2f mm, S=%.2f mm)", 螺距, scanner.参数.d, 层厚);
            
            if (参数已更改)
            {
                更新数据();
            }
        }
        
        // 3D绘图 - 使用ImPlot3D
        if (ImPlot3D::BeginPlot("螺旋轨迹与射线", ImVec2(-1, -1)))
        {
            // 绘制螺旋轨迹
            if (!螺旋轨迹点.empty())
            {
                std::vector<double> x_traj, y_traj, z_traj;
                x_traj.reserve(螺旋轨迹点.size());
                y_traj.reserve(螺旋轨迹点.size());
                z_traj.reserve(螺旋轨迹点.size());
                
                for (const auto& pt : 螺旋轨迹点)
                {
                    x_traj.push_back(pt.x);
                    y_traj.push_back(pt.y);
                    z_traj.push_back(pt.z);
                }
                
                ImPlot3D::PlotLine("螺旋轨迹", x_traj.data(), y_traj.data(), z_traj.data(), x_traj.size());
            }
            
            // 绘制当前角度的焦点位置
            Point3D focus = scanner.计算焦点位置(当前角度);
            double focus_x[1] = { focus.x };
            double focus_y[1] = { focus.y };
            double focus_z[1] = { focus.z };
            ImPlot3D::PlotScatter("焦点位置", focus_x, focus_y, focus_z, 1);
            
            // 绘制射线束（只绘制部分射线以避免性能问题）
            int ray_skip = std::max(1, static_cast<int>(当前射线束.size() / 50));  // 最多绘制50条射线
            int ray_count = 0;
            for (size_t i = 0; i < 当前射线束.size(); i += ray_skip)
            {
                const auto& ray = 当前射线束[i];
                double x_ray[2] = { ray.起点.x, ray.终点.x };
                double y_ray[2] = { ray.起点.y, ray.终点.y };
                double z_ray[2] = { ray.起点.z, ray.终点.z };
                
                char label[64];
                snprintf(label, sizeof(label), "射线_%d", ray_count++);
                ImPlot3D::PlotLine(label, x_ray, y_ray, z_ray, 2);
            }
            
            // 绘制旋转中心（随射线源z坐标移动）
            Point3D center = scanner.计算旋转中心(当前角度);
            double center_x[1] = { center.x };
            double center_y[1] = { center.y };
            double center_z[1] = { center.z };
            ImPlot3D::PlotScatter("旋转中心", center_x, center_y, center_z, 1);
            
            // 绘制切线平面
            if (显示切线平面 && !切线平面点.empty())
            {
                // 将平面点组织成网格进行绘制
                int resolution = 20;
                std::vector<double> plane_x, plane_y, plane_z;
                plane_x.reserve(切线平面点.size());
                plane_y.reserve(切线平面点.size());
                plane_z.reserve(切线平面点.size());
                
                for (const auto& pt : 切线平面点)
                {
                    plane_x.push_back(pt.x);
                    plane_y.push_back(pt.y);
                    plane_z.push_back(pt.z);
                }
                
                // 绘制平面网格线
                // 绘制行（沿切线方向）
                for (int i = 0; i < resolution; ++i)
                {
                    std::vector<double> row_x, row_y, row_z;
                    row_x.reserve(resolution);
                    row_y.reserve(resolution);
                    row_z.reserve(resolution);
                    
                    for (int j = 0; j < resolution; ++j)
                    {
                        int idx = i * resolution + j;
                        row_x.push_back(切线平面点[idx].x);
                        row_y.push_back(切线平面点[idx].y);
                        row_z.push_back(切线平面点[idx].z);
                    }
                    
                    char label[64];
                    snprintf(label, sizeof(label), "平面行_%d", i);
                    ImPlot3D::PlotLine(label, row_x.data(), row_y.data(), row_z.data(), resolution);
                }
                
                // 绘制列（沿径向方向）
                for (int j = 0; j < resolution; ++j)
                {
                    std::vector<double> col_x, col_y, col_z;
                    col_x.reserve(resolution);
                    col_y.reserve(resolution);
                    col_z.reserve(resolution);
                    
                    for (int i = 0; i < resolution; ++i)
                    {
                        int idx = i * resolution + j;
                        col_x.push_back(切线平面点[idx].x);
                        col_y.push_back(切线平面点[idx].y);
                        col_z.push_back(切线平面点[idx].z);
                    }
                    
                    char label[64];
                    snprintf(label, sizeof(label), "平面列_%d", j);
                    ImPlot3D::PlotLine(label, col_x.data(), col_y.data(), col_z.data(), resolution);
                }
            }
            
            ImPlot3D::EndPlot();
        }
        
        ImGui::End();
    }
};

int main()
{
    system("chcp 65001");

    // 创建螺旋锥束扫描器
    螺旋锥束扫描器 scanner;
    
    // 创建3D可视化渲染器
    debugger d(std::make_shared<spiral_3d_render>(scanner));
    return d.execute();
    
    // 如果需要使用原来的2D投影可视化，可以取消下面的注释
    /*
    锥形束正投影 projection;
    projection.display();
    projection.calculate();
    debugger d2(std::make_shared<plot_render>(projection));
    return d2.execute();
    */
}
