#include "Hybmesh.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <thread>
#include <mutex>
#include <vector>
#include <cmath>

// Функция уточнения [0, 1] -> [0, 1]
double reffun(double t) {
    return 1.0 + std::sin(0.5 * M_PI * (t - 1));
}

// Мьютекс для синхронизации вывода в консоль
std::mutex mtx;

void run_thread(int thread_id, const std::vector<double>& zcoords, int contour_points, int grid_size) {
    Hybmesh hm;

    // Создание контура
    auto contour = hm.add_circ_contour(Hybmesh::Point2(0.5, 0.5), 0.3, contour_points);

    // Вычисление координат блоков
    int row = thread_id / grid_size;
    int col = thread_id % grid_size;

    double x_start = static_cast<double>(col) / grid_size;
    double y_start = static_cast<double>(row) / grid_size;
    double x_end = static_cast<double>(col + 1) / grid_size;
    double y_end = static_cast<double>(row + 1) / grid_size;

    // Создание единичной квадратной сетки
    auto square = hm.add_unf_rect_grid({x_start, y_start}, {x_end, y_end}, 5, 5);

    // Объединение контура с квадратной сеткой
    auto res = hm.inscribe_grid(square, contour, "outside", 0);

    // Экструзия сетки
    auto extruded_grid = hm.extrude_grid(res, zcoords, 0, 0);

    // Экспорт сеточной модели
    std::string filename = "result_" + std::to_string(thread_id) + ".vtk";
    hm.export3d_grid_vtk(extruded_grid, filename);

    // Синхронизация доступа к консольному выводу
    std::lock_guard<std::mutex> lock(mtx);
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
        std::cerr << "Number of threads must be a perfect square (1, 4, 9, 16, ...)." << std::endl;
        return 1;
    }

    double minz = 0.0, maxz = 1.0;
    std::vector<double> zcoords(nz + 1);

    // Вычисление координат z с синусоидальным уточнением к z=0
    for (int i = 0; i <= nz; ++i) {
        double t = reffun(static_cast<double>(i) / nz);
        zcoords[i] = minz + (maxz - minz) * t;
    }

    // Измерение времени выполнения всех потоков
    auto start = std::chrono::high_resolution_clock::now();

    // Создание и запуск потоков
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.push_back(std::thread(run_thread, i, std::ref(zcoords), contour_points, grid_size));
    }

    // Ожидание завершения всех потоков
    for (auto& t : threads) {
        t.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_secs = end - start;
    std::cout << "Общее время выполнения: " << elapsed_secs.count() << " секунд" << std::endl;
    return 0;
}
