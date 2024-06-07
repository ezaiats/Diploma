#include "Hybmesh.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <omp.h>
#include <vector>

void process_grid(const Hybmesh::Point2& p1, const Hybmesh::Point2& p2, int thread_id, int contour_points, int nz) {
    Hybmesh hm;

    // Создание контура
    auto contour = hm.add_circ_contour(Hybmesh::Point2(0.5, 0.5), 0.3, contour_points);

    // Создание единичной квадратной сеткой
    auto square = hm.add_unf_rect_grid(p1, p2, 5, 5);

    // Объединение контура с единичной квадратной сеткой
    auto res = hm.inscribe_grid(square, contour, "outside", 0);

    // Вычисление координат z с синусоидальным уточнением к z=0
    double minz = 0.0, maxz = 1.0;
    std::vector<double> zcoords(nz + 1);
    for (int i = 0; i <= nz; ++i) {
        double t = static_cast<double>(i) / nz;
        t = 1.0 + std::sin(0.5 * M_PI * (t - 1));
        zcoords[i] = minz + (maxz - minz) * t;
    }

    // Экструзия сетки
    auto extruded_grid = hm.extrude_grid(res, zcoords, 0, 0);

    // Экспорт сеточной модели
    std::string filename = "result_" + std::to_string(thread_id) + ".vtk";
    hm.export3d_grid_vtk(extruded_grid, filename);

    std::cout << "Поток " << thread_id << " завершил работу." << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <number_of_threads> <contour_points> <nz>" << std::endl;
        return 1;
    }

    int num_threads = std::stoi(argv[1]);
    int contour_points = std::stoi(argv[2]);
    int nz = std::stoi(argv[3]);
    int grid_size = std::sqrt(num_threads);

    if (grid_size * grid_size != num_threads) {
        std::cerr << "Number of threads must be a perfect square (1, 4, 8, 16, ...)." << std::endl;
        return 1;
    }

    // Измерение времени выполнения всех потоков
    double start = omp_get_wtime();

    // Запуск потоков
    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        int row = thread_id / grid_size;
        int col = thread_id % grid_size;

        double x_start = static_cast<double>(col) / grid_size;
        double y_start = static_cast<double>(row) / grid_size;
        double x_end = static_cast<double>(col + 1) / grid_size;
        double y_end = static_cast<double>(row + 1) / grid_size;

        process_grid({x_start, y_start}, {x_end, y_end}, thread_id, contour_points, nz);
    }

    double end = omp_get_wtime();
    double elapsed_secs = end - start;

    std::cout << "Общее время выполнения: " << elapsed_secs << " секунд" << std::endl;
    return 0;
}
