#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <tuple>
#include <cmath>
#include <chrono>
#include <windows.h>

class Matrix {
public:
    struct SymbolicAction {
        int modified_row;
        int pivot_row;
        std::vector<int> affected_columns;  // Столбцы, которые будут изменены

        SymbolicAction(int i, int k, const std::vector<int>& cols)
            : modified_row(i), pivot_row(k), affected_columns(cols) {}
    };

    struct SymbolicPattern {
        std::vector<std::pair<int, int>> non_zero_positions;
        std::vector<int> unstable_columns;
        std::vector<SymbolicAction> actions;
    };

private:
    // LIL представление
    int rows;
    int cols;
    std::vector<std::vector<std::pair<int, long double>>> LIL;

    // CSR представление
    std::vector<long double> csr_data;
    std::vector<int> csr_col_indices;
    std::vector<int> csr_row_ptr;

    bool is_csr_current = false;
    SymbolicPattern symbolic_pattern;
    mutable std::vector<std::vector<long double>> last_dense;

    // Порог для определения численной неустойчивости
    const long double STABILITY_THRESHOLD = 1e-10;

    // Проверка на численную неустойчивость
    bool is_numerically_unstable(long double value) const {
        return std::abs(value) < STABILITY_THRESHOLD;
    }

    // Проверка индексов
    bool is_valid_index(int i, int j) const {
        return i >= 0 && i < rows && j >= 0 && j < cols;
    }

    // Обновить CSR представление
    void update_csr() {
        if (is_csr_current) return;

        csr_data.clear();
        csr_col_indices.clear();
        csr_row_ptr.resize(rows + 1, 0);

        for (int i = 0; i < rows; ++i) {
            for (const auto& p : LIL[i]) {
                if (p.first >= 0 && p.first < cols) {
                csr_data.push_back(p.second);
                csr_col_indices.push_back(p.first);
                csr_row_ptr[i + 1]++;
                }
            }
        }

        for (int i = 1; i <= rows; ++i) {
            csr_row_ptr[i] += csr_row_ptr[i - 1];
        }
        is_csr_current = true;
    }

    // CSR хелперы
    std::pair<int, int> csr_row_bounds(int row) const {
        if (row < 0 || row >= rows) {
            return { 0, 0 };
        }
        return { csr_row_ptr[row], csr_row_ptr[row + 1] };
    }

    long double* csr_find_element(int row, int col) {
        if (!is_valid_index(row, col)) return nullptr;
        
        auto start_end_pair = csr_row_bounds(row);
        auto it = std::lower_bound(csr_col_indices.begin() + start_end_pair.first,
            csr_col_indices.begin() + start_end_pair.second, col);
        return (it != csr_col_indices.end() && *it == col) ?
            &csr_data[it - csr_col_indices.begin()] : nullptr;
    }

public:
    Matrix(int r, int c) : rows(r), cols(c), csr_row_ptr(r + 1, 0) {
        if (r <= 0 || c <= 0) {
            throw std::invalid_argument("Matrix dimensions must be positive");
        }
    }

    // Сохранение в LIL
    void saveAsLIL(const std::vector<std::vector<long double>>& matrix) {
        if (matrix.size() != rows || matrix[0].size() != cols) {
            throw std::invalid_argument("Matrix dimensions do not match");
        }
        
        LIL.resize(rows);
        for (int i = 0; i < rows; ++i) {
            LIL[i].clear();
            for (int j = 0; j < cols; ++j) {
                if (matrix[i][j] != 0) {
                    LIL[i].emplace_back(j, matrix[i][j]);
                }
            }
        }
        is_csr_current = false;
    }

    // Сохранение в CSR
    void saveAsCSR(const std::vector<std::vector<long double>>& matrix) {
        saveAsLIL(matrix);
        update_csr();
    }

    // Конвертация в CSR
    void convertToCSR() {
        update_csr();
    }

    // Получение элемента
    long double get(int i, int j) const {
        if (!is_valid_index(i, j)) {
            throw std::out_of_range("Matrix index out of bounds");
        }

        if (is_csr_current && !csr_data.empty() && !csr_col_indices.empty() && !csr_row_ptr.empty()) {
            if (i >= csr_row_ptr.size() - 1) return 0;
            
            int start = csr_row_ptr[i];
            int end = csr_row_ptr[i + 1];
            
            if (start >= end || start >= csr_col_indices.size() || end > csr_col_indices.size()) {
                return 0;
            }

            auto it = std::lower_bound(csr_col_indices.begin() + start,
                csr_col_indices.begin() + end, j);
            
            if (it != csr_col_indices.end() && *it == j) {
                int index = it - csr_col_indices.begin();
                if (index < csr_data.size()) {
                    return csr_data[index];
                }
            }
            return 0;
        }

        for (const auto& p : LIL[i]) {
            if (p.first == j) return p.second;
        }
        return 0;
    }

    // Установка элемента (только для LIL)
    void set(int i, int j, long double value) {
        if (!is_valid_index(i, j)) {
            throw std::out_of_range("Matrix index out of bounds");
        }

        is_csr_current = false;
        
        // Убедимся, что LIL имеет нужный размер
        if (LIL.size() <= i) {
            LIL.resize(rows);
        }
        
        auto& row = LIL[i];
        for (auto& p : row) {
            if (p.first == j) {
                p.second = value;
                return;
            }
        }
        row.emplace_back(j, value);
        std::sort(row.begin(), row.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
    }

    // LU-разложение (объединенная версия)
    std::pair<Matrix, Matrix> luDecomposition(bool use_saved = false) {
        if (rows != cols) throw std::invalid_argument("Matrix must be square");

        if (!is_csr_current) {
            update_csr();
        }

        if (use_saved && !symbolic_pattern.actions.empty()) {
            return perform_fast_lu();
        }
        return perform_full_lu();
    }

    // Печать матрицы
    void printDense(int precision = 2) const {
        auto dense = getDense();
        for (const auto& row : dense) {
            for (long double val : row) {
                std::cout << std::setw(precision + 3)
                    << std::fixed << std::setprecision(precision)
                    << val << " ";
            }
            std::cout << "\n";
        }
    }

    // Умножение матриц
    Matrix multiply(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }
        
        Matrix result(rows, other.cols);
        
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                long double sum = 0.0;
                for (int k = 0; k < cols; ++k) {
                    sum += get(i, k) * other.get(k, j);
                }
                if (sum != 0) {
                    result.set(i, j, sum);
                }
            }
        }
        
        return result;
    }

    // Проверка LU-разложения
    bool checkLUDecomposition(const Matrix& L, const Matrix& U) const {
        Matrix product = L.multiply(U);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (std::abs(product.get(i, j) - get(i, j)) > 1e-10) {
                    return false;
                }
            }
        }
        return true;
    }

private:
    // Полное LU-разложение
    std::pair<Matrix, Matrix> perform_full_lu() {
        symbolic_pattern.actions.clear();
        symbolic_pattern.unstable_columns.clear();
        Matrix L(rows, cols);
        Matrix U = *this;

        // Устанавливаем диагональные элементы L равными 1
        for (int i = 0; i < rows; ++i) {
            L.set(i, i, 1.0);
        }

        for (int k = 0; k < rows; ++k) {
            long double pivot = U.get(k, k);
            
            // Проверка на неустойчивость
            if (is_numerically_unstable(pivot)) {
                symbolic_pattern.unstable_columns.push_back(k);
                continue;
            }

            // Собираем столбцы, которые будут изменены в текущем шаге
            std::vector<int> affected_cols;
            for (int j = k; j < cols; ++j) {
                if (U.get(k, j) != 0) {
                    affected_cols.push_back(j);
                }
            }

            // Классическое LU-разложение по Гауссу
            for (int i = k + 1; i < rows; ++i) {
                long double u_ik = U.get(i, k);
                if (std::abs(u_ik) < STABILITY_THRESHOLD) continue;

                // Вычисляем множитель
                long double multiplier = u_ik / pivot;
                L.set(i, k, multiplier);

                // Обновляем элементы в строке i
                for (int j : affected_cols) {
                    long double u_kj = U.get(k, j);
                    if (std::abs(u_kj) < STABILITY_THRESHOLD) continue;
                    
                    long double u_ij = U.get(i, j);
                    long double delta = multiplier * u_kj;
                    long double new_value = u_ij - delta;
                    
                    // Проверяем на численную устойчивость
                    if (std::abs(new_value) < STABILITY_THRESHOLD) {
                        U.set(i, j, 0.0);
                    } else {
                        U.set(i, j, new_value);
                    }
                }

                // Сохраняем действие для быстрого разложения
                symbolic_pattern.actions.emplace_back(i, k, affected_cols);
            }
        }

        return { L, U };
    }

    // Быстрое LU-разложение
    std::pair<Matrix, Matrix> perform_fast_lu() {
        Matrix L(rows, cols);
        Matrix U = *this;

        // Устанавливаем диагональные элементы L равными 1
        for (int i = 0; i < rows; ++i) {
            L.set(i, i, 1.0);
        }

        // Применяем все действия из символьной схемы
        for (const auto& action : symbolic_pattern.actions) {
            long double u_ik = U.get(action.modified_row, action.pivot_row);
            if (std::abs(u_ik) < STABILITY_THRESHOLD) continue;

            long double pivot = U.get(action.pivot_row, action.pivot_row);
            if (std::abs(pivot) < STABILITY_THRESHOLD) continue;

            long double multiplier = u_ik / pivot;
            L.set(action.modified_row, action.pivot_row, multiplier);

            // Обновляем только те элементы, которые указаны в affected_columns
            for (int j : action.affected_columns) {
                long double u_kj = U.get(action.pivot_row, j);
                if (std::abs(u_kj) < STABILITY_THRESHOLD) continue;
                
                long double u_ij = U.get(action.modified_row, j);
                long double delta = multiplier * u_kj;
                long double new_value = u_ij - delta;
                
                // Проверяем на численную устойчивость
                if (std::abs(new_value) < STABILITY_THRESHOLD) {
                    U.set(action.modified_row, j, 0.0);
                } else {
                    U.set(action.modified_row, j, new_value);
                }
            }
        }

        return { L, U };
    }

    // Получение плотного представления
    std::vector<std::vector<long double>> getDense() const {
        std::vector<std::vector<long double>> dense(
            rows, std::vector<long double>(cols, 0));

        if (is_csr_current && !csr_data.empty() && !csr_col_indices.empty() && !csr_row_ptr.empty()) {
            for (int i = 0; i < rows; ++i) {
                if (i >= csr_row_ptr.size() - 1) break;
                
                int start = csr_row_ptr[i];
                int end = csr_row_ptr[i + 1];
                
                if (start >= end || start >= csr_col_indices.size() || end > csr_col_indices.size()) {
                    continue;
                }

                for (int j = start; j < end; ++j) {
                    if (j >= csr_col_indices.size() || j >= csr_data.size()) break;
                    
                    int col = csr_col_indices[j];
                    if (col >= 0 && col < cols) {
                        dense[i][col] = csr_data[j];
                    }
                }
            }
        }
        else {
            for (int i = 0; i < rows && i < LIL.size(); ++i) {
                for (const auto& p : LIL[i]) {
                    if (p.first >= 0 && p.first < cols) {
                    dense[i][p.first] = p.second;
                    }
                }
            }
        }
        return dense;
    }
};

int main() {
    // Установка русского языка в консоль
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    setlocale(LC_ALL, "Russian");

    // Тест 1: Разреженная матрица с диагональным преобладанием
    Matrix A(10, 10);
    A.saveAsCSR({
        {2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0}
    });

    std::cout << "Тест 1: Разреженная матрица с диагональным преобладанием\n";
    std::cout << "Исходная матрица A:\n";
    A.printDense(6);

    // Полное разложение
    auto start_time = std::chrono::high_resolution_clock::now();
    auto ans = A.luDecomposition();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "\nВремя полного разложения: " << duration.count() << " микросекунд\n";
    
    // Быстрое разложение
    start_time = std::chrono::high_resolution_clock::now();
    auto ans_fast = A.luDecomposition(true);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Время быстрого разложения: " << duration.count() << " микросекунд\n";
    
    std::cout << "\nМатрица L:\n";
    ans.first.printDense(6);
    std::cout << "\nМатрица U:\n";
    ans.second.printDense(6);

    // Проверка LU-разложения
    std::cout << "\nПроверка LU-разложения (матрица A): " 
              << (A.checkLUDecomposition(ans.first, ans.second) ? "Успешно" : "Ошибка") << "\n";
    
    // После быстрого разложения матрицы A:
    std::cout << "Проверка быстрого LU-разложения (матрица A): " 
              << (A.checkLUDecomposition(ans_fast.first, ans_fast.second) ? "Успешно" : "Ошибка") << "\n";

    // Тест 2: Разреженная матрица с несколькими ненулевыми элементами
    Matrix B(10, 10);
    B.saveAsCSR({
        {2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0}
    });

    std::cout << "\nТест 2: Разреженная матрица с несколькими ненулевыми элементами\n";
    std::cout << "Исходная матрица B:\n";
    B.printDense(6);

    // Полное разложение
    start_time = std::chrono::high_resolution_clock::now();
    auto ans2 = B.luDecomposition();
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "\nВремя полного разложения: " << duration.count() << " микросекунд\n";
    
    // Быстрое разложение
    start_time = std::chrono::high_resolution_clock::now();
    auto ans2_fast = B.luDecomposition(true);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Время быстрого разложения: " << duration.count() << " микросекунд\n";
    
    std::cout << "\nМатрица L:\n";
    ans2.first.printDense(6);
    std::cout << "\nМатрица U:\n";
    ans2.second.printDense(6);

    // Проверка LU-разложения
    std::cout << "\nПроверка LU-разложения (матрица B): " 
              << (B.checkLUDecomposition(ans2.first, ans2.second) ? "Успешно" : "Ошибка") << "\n";
    
    // После быстрого разложения матрицы B:
    std::cout << "Проверка быстрого LU-разложения (матрица B): " 
              << (B.checkLUDecomposition(ans2_fast.first, ans2_fast.second) ? "Успешно" : "Ошибка") << "\n";

    // Тест 3: Разреженная матрица с численно неустойчивыми элементами
    Matrix C(10, 10);
    C.saveAsCSR({
        {2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 2.0, 1e-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0}
    });

    std::cout << "\nТест 3: Разреженная матрица с численно неустойчивыми элементами\n";
    std::cout << "Исходная матрица C:\n";
    C.printDense(6);

    // Полное разложение
    start_time = std::chrono::high_resolution_clock::now();
    auto ans3 = C.luDecomposition();
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "\nВремя полного разложения: " << duration.count() << " микросекунд\n";
    
    // Быстрое разложение
    start_time = std::chrono::high_resolution_clock::now();
    auto ans3_fast = C.luDecomposition(true);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Время быстрого разложения: " << duration.count() << " микросекунд\n";
    
    std::cout << "\nМатрица L:\n";
    ans3.first.printDense(6);
    std::cout << "\nМатрица U:\n";
    ans3.second.printDense(6);

    // Проверка LU-разложения
    std::cout << "\nПроверка LU-разложения (матрица C): " 
              << (C.checkLUDecomposition(ans3.first, ans3.second) ? "Успешно" : "Ошибка") << "\n";
    
    // После быстрого разложения матрицы C:
    std::cout << "Проверка быстрого LU-разложения (матрица C): " 
              << (C.checkLUDecomposition(ans3_fast.first, ans3_fast.second) ? "Успешно" : "Ошибка") << "\n";

    // Тест 4: Более заполненная разреженная матрица
    Matrix D(10, 10);
    D.saveAsCSR({
        {2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0}
    });

    std::cout << "\nТест 4: Более заполненная разреженная матрица\n";
    std::cout << "Исходная матрица D:\n";
    D.printDense(6);

    // Полное разложение
    start_time = std::chrono::high_resolution_clock::now();
    auto ans4 = D.luDecomposition();
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "\nВремя полного разложения: " << duration.count() << " микросекунд\n";
    
    std::cout << "\nМатрица L (полное разложение):\n";
    ans4.first.printDense(6);
    std::cout << "\nМатрица U (полное разложение):\n";
    ans4.second.printDense(6);
    
    // Быстрое разложение
    start_time = std::chrono::high_resolution_clock::now();
    auto ans4_fast = D.luDecomposition(true);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "\nВремя быстрого разложения: " << duration.count() << " микросекунд\n";
    
    std::cout << "\nМатрица L (быстрое разложение):\n";
    ans4_fast.first.printDense(6);
    std::cout << "\nМатрица U (быстрое разложение):\n";
    ans4_fast.second.printDense(6);

    // Проверка LU-разложения
    std::cout << "\nПроверка LU-разложения (матрица D): " 
              << (D.checkLUDecomposition(ans4.first, ans4.second) ? "Успешно" : "Ошибка") << "\n";
    
    // После быстрого разложения матрицы D:
    std::cout << "Проверка быстрого LU-разложения (матрица D): " 
              << (D.checkLUDecomposition(ans4_fast.first, ans4_fast.second) ? "Успешно" : "Ошибка") << "\n";

    // Тест 5: Сравнение производительности для матриц с одинаковой структурой
    std::cout << "\nТест 5: Сравнение производительности для матриц с одинаковой структурой\n";
    
    // Первая матрица (как в тесте 4)
    Matrix E(10, 10);
    E.saveAsCSR({
        {2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0}
    });

    std::cout << "Первая матрица E:\n";
    E.printDense(6);

    // Полное разложение первой матрицы
    start_time = std::chrono::high_resolution_clock::now();
    auto ans5 = E.luDecomposition();
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "\nВремя полного разложения первой матрицы: " << duration.count() << " микросекунд\n";
    
    std::cout << "\nМатрица L (первая матрица, полное разложение):\n";
    ans5.first.printDense(6);
    std::cout << "\nМатрица U (первая матрица, полное разложение):\n";
    ans5.second.printDense(6);

    // Проверка LU-разложения
    std::cout << "\nПроверка LU-разложения (матрица E): " 
              << (E.checkLUDecomposition(ans5.first, ans5.second) ? "Успешно" : "Ошибка") << "\n";

    // Вторая матрица с той же структурой, но другими данными
    Matrix F(10, 10);
    F.saveAsCSR({
        {3.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {2.0, 3.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 2.0, 3.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0},
        {2.0, 0.0, 2.0, 3.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0},
        {0.0, 2.0, 0.0, 2.0, 3.0, 2.0, 0.0, 2.0, 0.0, 0.0},
        {0.0, 0.0, 2.0, 0.0, 2.0, 3.0, 2.0, 0.0, 2.0, 0.0},
        {0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 3.0, 2.0, 0.0, 2.0},
        {0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 3.0, 2.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 3.0, 2.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 3.0}
    });

    std::cout << "\nВторая матрица F:\n";
    F.printDense(6);

    // Быстрое разложение второй матрицы, используя структуру от матрицы E
    start_time = std::chrono::high_resolution_clock::now();
    auto ans6_fast = F.luDecomposition(true);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "\nВремя быстрого разложения второй матрицы: " << duration.count() << " микросекунд\n";
    
    std::cout << "\nМатрица L (вторая матрица, быстрое разложение):\n";
    ans6_fast.first.printDense(6);
    std::cout << "\nМатрица U (вторая матрица, быстрое разложение):\n";
    ans6_fast.second.printDense(6);

    // Проверка быстрого LU-разложения
    std::cout << "Проверка быстрого LU-разложения (матрица F): " 
              << (F.checkLUDecomposition(ans6_fast.first, ans6_fast.second) ? "Успешно" : "Ошибка") << "\n";

    // Полное разложение второй матрицы (для сравнения)
    start_time = std::chrono::high_resolution_clock::now();
    auto ans6 = F.luDecomposition();
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "\nВремя полного разложения второй матрицы: " << duration.count() << " микросекунд\n";
    
    std::cout << "\nМатрица L (вторая матрица, полное разложение):\n";
    ans6.first.printDense(6);
    std::cout << "\nМатрица U (вторая матрица, полное разложение):\n";
    ans6.second.printDense(6);

    // Проверка полного LU-разложения
    std::cout << "\nПроверка полного LU-разложения (матрица F): " 
              << (F.checkLUDecomposition(ans6.first, ans6.second) ? "Успешно" : "Ошибка") << "\n";

    return 0;
}