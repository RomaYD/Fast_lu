#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include <map>
#include <windows.h>
#include <fstream>
#include <chrono>
#include <limits>
#include <unordered_map>
#include <climits>


// Структура для хранения информации о численной неустойчивости
struct struct_of_decomp {
    std::vector<std::vector<long long>> data;
    std::vector<long long> csr_col;
    std::vector<long long> csr_row;
};

struct NumericalInstability {
    int k;              // Номер шага
    int i;              // Номер строки
    int j;              // Номер столбца
    long double value;  // Значение, вызвавшее неустойчивость
    long double threshold; // Порог устойчивости

    NumericalInstability(int _k, int _i, int _j, long double _value, long double _threshold)
        : k(_k), i(_i), j(_j), value(_value), threshold(_threshold) {}
};

class Matrix {
private:
    // Размеры матрицы
    int rows;
    int cols;

    // LIL представление - список строк, каждая строка - список пар (индекс, значение)
    std::vector<std::list<std::pair<int, long double>>> LIL;

    // CSR представление
    std::vector<long double> csr_data;      // Значения ненулевых элементов
    std::vector<long long> csr_col_indices;       // Индексы столбцов
    std::vector<long long> csr_row_ptr;           // Указатели на начало строк

    // Проверка индексов
    bool is_valid_index(int i, int j) const {
        return i >= 0 && i < rows && j >= 0 && j < cols;
    }

public:
    // Конструктор
    Matrix(int r, int c) : rows(r), cols(c) {
        if (r <= 0 || c <= 0) {
            throw std::invalid_argument("Размеры матрицы должны быть положительными");
        }
        LIL.resize(rows);
    }

    // Получение элемента
    long double get(int i, int j) const {
        if (!is_valid_index(i, j)) {
            throw std::out_of_range("Индекс матрицы вне границ");
        }

        // Простой поиск по списку, как в статье
        for (const auto& p : LIL[i]) {
            if (p.first == j) return p.second;
        }
        return 0;
    }

    // Установка элемента
    void set(int i, int j, long double value) {
        if (!is_valid_index(i, j)) {
            throw std::out_of_range("Индекс матрицы вне границ");
        }

        auto& row = LIL[i];
        
        // Поиск позиции для вставки или обновления
        auto it = row.begin();
        while (it != row.end()) {
            if (it->first > j) break;
            else if (it->first == j) {
                if (std::abs(value) > 1e-10) {
                    it->second = value;
                } else {
                    row.erase(it);
                }
                return;
            }
            ++it;
        }

        // Вставка нового элемента
        if (std::abs(value) > 1e-10) {
            row.insert(it, std::make_pair(j, value));
        }
    }

    // Сохранение матрицы из плотного формата
    void saveFromDense(const std::vector<std::vector<long double>>& matrix) {
        if (matrix.size() != rows || matrix[0].size() != cols) {
            throw std::invalid_argument("Размеры матрицы не совпадают");
        }
        
        LIL.resize(rows);
        for (int i = 0; i < rows; ++i) {
            LIL[i].clear();
            for (int j = 0; j < cols; ++j) {
                if (std::abs(matrix[i][j]) > 1e-10) {
                    LIL[i].push_back(std::make_pair(j, matrix[i][j]));
                }
            }
        }
    }

    // Перевод из LIL в плотное представление
    std::vector<std::vector<long double>> toDense() const {
        std::vector<std::vector<long double>> dense(rows, std::vector<long double>(cols, 0.0));
        for (int i = 0; i < rows; ++i) {
            for (const auto& p : LIL[i]) {
                dense[i][p.first] = p.second;
            }
        }
        return dense;
    }

    // Перевод из CSR в плотное представление
    std::vector<std::vector<long double>> toDenseFromCSR() const {
        std::vector<std::vector<long double>> dense(rows, std::vector<long double>(cols, 0.0));
        for (int i = 0; i < rows; ++i) {
            int start = csr_row_ptr[i];
            int end = csr_row_ptr[i + 1];
            for (int k = start; k < end; ++k) {
                dense[i][csr_col_indices[k]] = csr_data[k];
            }
        }
        return dense;
    }

    // Прямой перевод из LIL в CSR
    void toCSR() {
        csr_data.clear();
        csr_col_indices.clear();
        csr_row_ptr.clear();
        csr_row_ptr.push_back(0);  // Начало первой строки

        for (int i = 0; i < rows; ++i) {
            for (const auto& p : LIL[i]) {
                csr_data.push_back(p.second);
                csr_col_indices.push_back(p.first);
            }
            csr_row_ptr.push_back(csr_data.size());  // Начало следующей строки
        }
    }

    // Сохранение в CSR формат
    void saveAsCSR() {
        toCSR();
    }

    // Получение элемента из CSR формата
    long double getFromCSR(int i, int j) const {
        if (!is_valid_index(i, j)) {
            throw std::out_of_range("Индекс матрицы вне границ");
        }

        // Ищем элемент в строке i
        int start = csr_row_ptr[i];
        int end = csr_row_ptr[i + 1];
        
        for (int k = start; k < end; ++k) {
            if (csr_col_indices[k] == j) {
                return csr_data[k];
            }
        }
        return 0;
    }

    // Печать матрицы в плотном формате
    void printDense(int precision = 2) const {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                long double val = get(i, j);
                if (std::abs(val) < 1e-10) val = 0;
                std::cout << std::setw(precision + 3)
                    << std::fixed << std::setprecision(precision)
                    << val << " ";
            }
            std::cout << "\n";
        }
    }

    // Печать матрицы в формате LIL
    void printLIL() const {
        std::cout << "LIL формат:\n";
        for (int i = 0; i < rows; ++i) {
            std::cout << "Строка " << i << ": ";
            for (const auto& p : LIL[i]) {
                std::cout << "(" << p.first << ", " << p.second << ") ";
            }
            std::cout << "\n";
        }
    }

    // Печать матрицы в формате CSR
    void printCSR() const {
        std::cout << "CSR формат:\n";
        std::cout << "data: ";
        for (const auto& val : csr_data) {
            std::cout << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\ncol_indices: ";
        for (const auto& idx : csr_col_indices) {
            std::cout << idx << " ";
        }
        std::cout << "\nrow_ptr: ";
        for (const auto& ptr : csr_row_ptr) {
            std::cout << ptr << " ";
        }
        std::cout << "\n";
    }

    // Получение размера матрицы
    std::pair<int, int> getSize() const {
        return { rows, cols };
    }

    // LU разложение с сохранением шагов
    friend std::pair<std::pair<Matrix, Matrix>, std::pair<struct_of_decomp, struct_of_decomp>> LUDecomposition(const Matrix& matrix) {
        if (matrix.rows != matrix.cols) {
            throw std::invalid_argument("Матрица должна быть квадратной для LU разложения");
        }

        Matrix L(matrix.rows, matrix.cols);
        Matrix U(matrix.rows, matrix.cols);
        std::pair<struct_of_decomp, struct_of_decomp> struc;
        std::vector<std::vector<long long>> struct_of_decomposition_L(3);
        std::vector<std::vector<long long>> struct_of_decomposition_U(3);

        // Инициализация L единичной матрицей
        for (int i = 0; i < matrix.rows; ++i) {
            L.csr_data.push_back(1.0);
            L.csr_col_indices.push_back(i);
            L.csr_row_ptr.push_back(i);
        }
        L.csr_row_ptr.push_back(matrix.rows);

        // Копируем исходную матрицу в U
        U = matrix;  // Используем оператор присваивания для копирования всех данных
        long long pivot_ind_U;
        // LU разложение
        for (int k = 0; k < matrix.rows - 1; ++k) {
            // Получаем опорный элемент U[k][k]
            long double pivot = 0;
            for (int j = U.csr_row_ptr[k]; j < U.csr_row_ptr[k + 1]; ++j) {
                if (U.csr_col_indices[j] == k) {
                    pivot = U.csr_data[j];
                    pivot_ind_U = j;
                    break;
                }
            }

            // Если опорный элемент слишком мал, используем небольшое значение
            
            if (std::abs(pivot) < 1e-10) {
                pivot = 1e-10;
                // Обновляем pivot в U
                bool found = false;
                for (int j = U.csr_row_ptr[k]; j < U.csr_row_ptr[k + 1]; ++j) {
                    if (U.csr_col_indices[j] == k) {
                        U.csr_data[j] = pivot;
                        pivot_ind_U = j;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    // Добавляем новый элемент в U
                    int insert_pos = U.csr_row_ptr[k];
                    while (insert_pos < U.csr_row_ptr[k + 1] && U.csr_col_indices[insert_pos] < k) {
                        insert_pos++;
                    }
                    struc.second.data.push_back({ insert_pos, pivot_ind_U });
                    U.csr_data.insert(U.csr_data.begin() + insert_pos, pivot);
                    U.csr_col_indices.insert(U.csr_col_indices.begin() + insert_pos, k);
                    for (int row = k + 1; row <= matrix.rows; ++row) {
                        U.csr_row_ptr[row]++;
                    }
                }
            }
            int j_of_u_kj;
            int j_of_u_ik;
            // Обрабатываем каждую строку ниже k
            for (int i = k + 1; i < matrix.rows; ++i) {
                // Получаем элемент U[i][k]
                long double u_ik = 0;
                for (int j = U.csr_row_ptr[i]; j < U.csr_row_ptr[i + 1]; ++j) {
                    if (U.csr_col_indices[j] == k) {
                        u_ik = U.csr_data[j];
                        j_of_u_ik = j;
                        break;
                    }
                }

                if (std::abs(u_ik) > 1e-10) {
                    long double factor = u_ik / pivot;

                    // Добавляем элемент в L
                    int insert_pos = L.csr_row_ptr[i];
                    while (insert_pos < L.csr_row_ptr[i + 1] && L.csr_col_indices[insert_pos] < k) {
                        insert_pos++;
                    }
                    struc.first.data.push_back({insert_pos, pivot_ind_U, j_of_u_ik });
                    L.csr_data.insert(L.csr_data.begin() + insert_pos, factor);
                    L.csr_col_indices.insert(L.csr_col_indices.begin() + insert_pos, k);
                    for (int row = i + 1; row <= matrix.rows; ++row) {
                        L.csr_row_ptr[row]++;
                    }

                    // Обновляем строку i в матрице U
                    for (int j = U.csr_row_ptr[k]; j < U.csr_row_ptr[k + 1]; ++j) {
                        int col = U.csr_col_indices[j];
                        if (col >= k) {
                            long double u_kj = U.csr_data[j];
                            j_of_u_kj = j;
                            
                            // Ищем или создаем элемент U[i][col]
                            bool found = false;
                            for (int m = U.csr_row_ptr[i]; m < U.csr_row_ptr[i + 1]; ++m) {
                                if (U.csr_col_indices[m] == col) {
                                    U.csr_data[m] -= factor * u_kj;
                                    struc.second.data.push_back({m, pivot_ind_U, j_of_u_ik, j_of_u_kj, 1});
                                    found = true;
                                    break;
                                }
                            }
                            
                            if (!found && std::abs(factor * u_kj) > 1e-10) {
                                // Добавляем новый элемент
                                int insert_pos = U.csr_row_ptr[i];
                                while (insert_pos < U.csr_row_ptr[i + 1] && U.csr_col_indices[insert_pos] < col) {
                                    insert_pos++;
                                }
                                struc.second.data.push_back({insert_pos, pivot_ind_U, j_of_u_ik, j_of_u_kj, 0});
                                U.csr_data.insert(U.csr_data.begin() + insert_pos, -factor * u_kj);
                                U.csr_col_indices.insert(U.csr_col_indices.begin() + insert_pos, col);
                                for (int row = i + 1; row <= matrix.rows; ++row) {
                                    U.csr_row_ptr[row]++;
                                }
                            }
                        }
                        if (struc.second.data.size() == 6)
                            int a = 1 + 1;
                    }

                    // Удаляем элемент U[i][k]
                    for (int j = U.csr_row_ptr[i]; j < U.csr_row_ptr[i + 1]; ++j) {
                        if (U.csr_col_indices[j] == k) {
                            struc.second.data.push_back({j});
                            U.csr_data.erase(U.csr_data.begin() + j);
                            U.csr_col_indices.erase(U.csr_col_indices.begin() + j);
                            for (int row = i + 1; row <= matrix.rows; ++row) {
                                U.csr_row_ptr[row]--;
                            }
                            break;
                        }
                    }
                }
            }
        }
        struct_of_decomposition_L[1] = L.csr_col_indices;
        struct_of_decomposition_U[1] = U.csr_col_indices;
        struct_of_decomposition_L[2] = L.csr_row_ptr;
        struct_of_decomposition_U[2] = U.csr_row_ptr;
        
        struc.first.csr_col = L.csr_col_indices;
        struc.first.csr_row = L.csr_row_ptr;
        struc.second.csr_col = U.csr_col_indices;
        struc.second.csr_row = U.csr_row_ptr;
        return {{L, U}, struc};
    }

    // Быстрое LU разложение с использованием структуры предыдущего разложения
    friend std::pair<Matrix, Matrix> FastLUDecomposition(const Matrix& matrix, const std::pair<struct_of_decomp, struct_of_decomp> previous_steps) 
    {
        Matrix L(matrix.rows, matrix.cols);
        Matrix U(matrix.rows, matrix.cols);
        long double u_ik;
        int prev_ind_u_ik = INT_MAX;
        long double L_insert = 0;
        long double factor;
        long double pivot;
        for (int i = 0; i < matrix.rows; ++i) {
            L.csr_data.push_back(1.0);
            L.csr_col_indices.push_back(i);
            L.csr_row_ptr.push_back(i);
        }
        L.csr_row_ptr.push_back(matrix.rows);
        U.csr_data = matrix.csr_data;
        L.csr_row_ptr = previous_steps.first.csr_row;
        L.csr_col_indices = previous_steps.first.csr_col;
        U.csr_col_indices = previous_steps.second.csr_col;
        U.csr_row_ptr = previous_steps.second.csr_row;
        // Идём по структуре и вычисляем элементы
        // todo написать более подробные коменты
        for (int i = 0; i < previous_steps.second.data.size(); i++) {
            if (i == 25)
                int a = 1 + 1;
            if (previous_steps.second.data[i].size() == 1) {
                U.csr_data.erase(U.csr_data.begin() + static_cast<int>(previous_steps.second.data[i][0]));
                prev_ind_u_ik = INT_MAX;
            }
            else
            {
                if (previous_steps.second.data[i].size() == 2) {
                    pivot = U.csr_data[static_cast<int>(previous_steps.second.data[i][1])];
                    U.csr_data.insert(U.csr_data.begin() + static_cast<int>(previous_steps.second.data[i][0]), pivot);
                }
                else
                {
                    if (static_cast<int>(previous_steps.second.data[i][4] == 1)) {
                        if (prev_ind_u_ik == INT_MAX) {
                            prev_ind_u_ik = static_cast<int>(previous_steps.second.data[i][2]);
                            u_ik = U.csr_data[prev_ind_u_ik];
                            pivot = U.csr_data[static_cast<int>(previous_steps.second.data[i][1])];
                            factor = u_ik / pivot;
                            L.csr_data.insert(L.csr_data.begin() + static_cast<int>(previous_steps.first.data[L_insert][0]), factor);
                            L_insert++;
                        }
                        if (prev_ind_u_ik != static_cast<int>(previous_steps.second.data[i][2])) {
                            prev_ind_u_ik = static_cast<int>(previous_steps.second.data[i][2]);
                            u_ik = U.csr_data[prev_ind_u_ik];
                            pivot = U.csr_data[static_cast<int>(previous_steps.second.data[i][1])];
                            factor = u_ik / pivot;
                            L.csr_data.insert(L.csr_data.begin() + static_cast<int>(previous_steps.first.data[L_insert][0]), factor);
                            L_insert++;
                        }
                        U.csr_data[static_cast<int>(previous_steps.second.data[i][0])] -= factor * U.csr_data[static_cast<int>(previous_steps.second.data[i][3])];
                    }
                    else
                    {
                        if (prev_ind_u_ik == INT_MAX) {
                            prev_ind_u_ik = static_cast<int>(previous_steps.second.data[i][2]);
                            u_ik = U.csr_data[prev_ind_u_ik];
                            pivot = U.csr_data[static_cast<int>(previous_steps.second.data[i][1])];
                            factor = u_ik / pivot;
                            L.csr_data.insert(L.csr_data.begin() + static_cast<int>(previous_steps.first.data[L_insert][0]), factor);
                            L_insert++;
                        }
                        if (prev_ind_u_ik != static_cast<int>(previous_steps.second.data[i][2])) {
                            prev_ind_u_ik = static_cast<int>(previous_steps.second.data[i][2]);
                            u_ik = U.csr_data[prev_ind_u_ik];
                            pivot = U.csr_data[static_cast<int>(previous_steps.second.data[i][1])];
                            factor = u_ik / pivot;
                            L.csr_data.insert(L.csr_data.begin() + static_cast<int>(previous_steps.first.data[L_insert][0]), factor);
                            L_insert++;
                        }
                        U.csr_data.insert(U.csr_data.begin() + static_cast<int>(previous_steps.second.data[i][0]), -factor * U.csr_data[static_cast<int>(previous_steps.second.data[i][3])]);
                        
                    }
                }
            }
        }
        return { L, U };

    }

    // Умножение матриц
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Несовместимые размеры матриц для умножения");
        }

        Matrix result(rows, other.cols);
        
            for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                long double sum = 0;
                for (int k = 0; k < cols; ++k) {
                    sum += get(i, k) * other.get(k, j);
                }
                if (std::abs(sum) > 1e-10) {
                    result.set(i, j, sum);
                }
            }
        }
        
        return result;
    }

    // Проверка равенства матриц
    bool operator==(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            return false;
        }

            for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (std::abs(get(i, j) - other.get(i, j)) > 1e-10) {
                    return false;
                }
            }
        }
        return true;
    }

    // Перевод матрицы в LIL формат
    void toLIL() {
        LIL.clear();
        LIL.resize(rows);
        
        // Проходим по всем элементам в CSR формате
        for (int i = 0; i < rows; ++i) {
            int start = csr_row_ptr[i];
            int end = csr_row_ptr[i + 1];
            
            // Добавляем все ненулевые элементы из строки i
            for (int k = start; k < end; ++k) {
                int j = csr_col_indices[k];
                long double value = csr_data[k];
                if (std::abs(value) > 1e-10) {
                    LIL[i].push_back(std::make_pair(j, value));
                }
            }
            
            // Сортируем элементы по индексу столбца
            LIL[i].sort([](const auto& a, const auto& b) {
                return a.first < b.first;
            });
        }
    }
};

// Тестовая функция
void testMatrix() {
    std::cout << "Тест матрицы 5x5:\n\n";

    // Создаем разреженную матрицу 5x5
    Matrix A(5, 5);
    A.saveFromDense({
        {1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 2.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 3.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 4.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 5.0}
    });
    A.toCSR();
    std::cout << "Матрица A (диагональная):\n";
    A.printDense(6);
    std::cout << "\n";

    // LU разложение
    std::cout << "LU разложение матрицы A:\n";
    auto start_time = std::chrono::high_resolution_clock::now();
    std::pair<std::pair<Matrix, Matrix>, std::pair<struct_of_decomp, struct_of_decomp>> result = LUDecomposition(A);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Время выполнения обычного LU разложения: " << duration.count() << " микросекунд\n\n";
    
    Matrix& L = result.first.first;
    Matrix& U = result.first.second;
    
    std::cout << "Матрица L (плотный формат):\n";
    auto L_dense = L.toDenseFromCSR();
    for (const auto& row : L_dense) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    std::cout << "Матрица U (плотный формат):\n";
    auto U_dense = U.toDenseFromCSR();
    for (const auto& row : U_dense) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // Проверка L * U = A
    L.toLIL();
    U.toLIL();
    Matrix LU = L * U;
    std::cout << "Результат умножения L * U:\n";
    auto LU_dense = LU.toDense();
    for (const auto& row : LU_dense) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    A.toLIL();
    std::cout << "Матрицы A и L*U " << (A == LU ? "совпадают" : "не совпадают") << "\n\n";

    // Создаем разреженную матрицу 10x10
    Matrix B(10, 10);
    std::vector<std::vector<long double>> B_dense(10, std::vector<long double>(10, 0.0));
    
    // Заполняем матрицу B с определенной структурой
    for (int i = 0; i < 10; ++i) {
        // Диагональные элементы
        B_dense[i][i] = 1.0 + i % 5;  // Значения от 1 до 5
        
        // Элементы на расстоянии 2 от диагонали
        if (i + 2 < 10) B_dense[i][i + 2] = 0.5;
        if (i - 2 >= 0) B_dense[i][i - 2] = 0.5;
        
        // Элементы на расстоянии 3 от диагонали
        if (i + 3 < 10) B_dense[i][i + 3] = 0.3;
        if (i - 3 >= 0) B_dense[i][i - 3] = 0.3;
    }
    
    B.saveFromDense(B_dense);
    B.toCSR();

    std::cout << "Матрица B (разреженная 10x10):\n";
    B.printDense(6);
    std::cout << "\n";

    // LU разложение
    std::cout << "LU разложение матрицы B:\n";
    start_time = std::chrono::high_resolution_clock::now();
    std::pair<std::pair<Matrix, Matrix>, std::pair<struct_of_decomp, struct_of_decomp>> result_B = LUDecomposition(B);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Время выполнения обычного LU разложения: " << duration.count() << " микросекунд\n\n";
    
    Matrix& L_B = result_B.first.first;
    Matrix& U_B = result_B.first.second;
    
    std::cout << "Матрица L (плотный формат):\n";
    auto L_dense_B = L_B.toDenseFromCSR();
    for (const auto& row : L_dense_B) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    std::cout << "Матрица U (плотный формат):\n";
    auto U_dense_B = U_B.toDenseFromCSR();
    for (const auto& row : U_dense_B) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // Проверка L * U = B
    L_B.toLIL();
    U_B.toLIL();
    Matrix LU_B = L_B * U_B;
    std::cout << "Результат умножения L * U:\n";
    auto LU_dense_B = LU_B.toDense();
    for (const auto& row : LU_dense_B) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    B.toLIL();
    std::cout << "Матрицы B и L*U " << (B == LU_B ? "совпадают" : "не совпадают") << "\n";

    // Тест быстрого LU разложения
    std::cout << "\nТест быстрого LU разложения матрицы B:\n";
    start_time = std::chrono::high_resolution_clock::now();
    std::pair<Matrix, Matrix> fast_result = FastLUDecomposition(B, result_B.second);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Время выполнения быстрого LU разложения: " << duration.count() << " микросекунд\n\n";
    
    Matrix& L_fast = fast_result.first;
    Matrix& U_fast = fast_result.second;
    //const std::vector<NumericalInstability>& instabilities = fast_result.second;
    
    std::cout << "Матрица L (плотный формат):\n";
    auto L_dense_fast = L_fast.toDenseFromCSR();
    for (const auto& row : L_dense_fast) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    std::cout << "Матрица U (плотный формат):\n";
    auto U_dense_fast = U_fast.toDenseFromCSR();
    for (const auto& row : U_dense_fast) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    /*if (!instabilities.empty()) {
        std::cout << "Обнаружены места численной неустойчивости:\n";
        for (const auto& inst : instabilities) {
            std::cout << "Шаг " << inst.k << ": "
                      << "Строка " << inst.i << ", "
                      << "Столбец " << inst.j << ", "
                      << "Значение = " << inst.value << ", "
                      << "Порог = " << inst.threshold << "\n";
        }
        std::cout << "\n";
    }*/

    // Проверка L * U = B для быстрого разложения
    L_fast.toLIL();
    U_fast.toLIL();
    Matrix LU_fast = L_fast * U_fast;
    std::cout << "Результат умножения L * U (быстрое разложение):\n";
    auto LU_dense_fast = LU_fast.toDense();
    for (const auto& row : LU_dense_fast) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "Матрицы B и L*U (быстрое разложение) " 
              << (B == LU_fast ? "совпадают" : "не совпадают") << "\n\n";

    // Создаем матрицу C с такой же структурой как у B
    Matrix C(10, 10);
    std::vector<std::vector<long double>> C_dense(10, std::vector<long double>(10, 0.0));
    
    // Заполняем матрицу C с той же структурой, но другими значениями
    for (int i = 0; i < 10; ++i) {
        // Диагональные элементы
        C_dense[i][i] = 2.0 + i % 5;  // Значения от 2 до 6
        
        // Элементы на расстоянии 2 от диагонали
        if (i + 2 < 10) C_dense[i][i + 2] = 0.7;
        if (i - 2 >= 0) C_dense[i][i - 2] = 0.7;
        
        // Элементы на расстоянии 3 от диагонали
        if (i + 3 < 10) C_dense[i][i + 3] = 0.4;
        if (i - 3 >= 0) C_dense[i][i - 3] = 0.4;
    }
    
    C.saveFromDense(C_dense);
    C.toCSR();
    std::cout << "Матрица C (та же структура что и у B):\n";
    C.printDense(6);
    std::cout << "\n";

    // Быстрое LU разложение матрицы C с использованием структуры B
    std::cout << "Быстрое LU разложение матрицы C (используя структуру B):\n";
    start_time = std::chrono::high_resolution_clock::now();
    std::pair<Matrix, Matrix> fast_result_C = FastLUDecomposition(C, result_B.second);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Время выполнения быстрого LU разложения: " << duration.count() << " микросекунд\n\n";
    
    Matrix& L_fast_C = fast_result_C.first;
    Matrix& U_fast_C = fast_result_C.second;
    //const std::vector<NumericalInstability>& instabilities_C = fast_result_C.second;

    // Обычное LU разложение матрицы C
    std::cout << "Обычное LU разложение матрицы C:\n";
    start_time = std::chrono::high_resolution_clock::now();
    std::pair<std::pair<Matrix, Matrix>, std::pair<struct_of_decomp, struct_of_decomp>> result_C = LUDecomposition(C);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "Время выполнения обычного LU разложения: " << duration.count() << " микросекунд\n\n";
    
    Matrix& L_C = result_C.first.first;
    Matrix& U_C = result_C.first.second;

    // Проверка результатов
    L_C.toLIL();
    U_C.toLIL();
    Matrix LU_C = L_C * U_C;
    L_fast_C.toLIL();
    U_fast_C.toLIL();
    Matrix LU_fast_C = L_fast_C * U_fast_C;

    std::cout << "Результат умножения L * U (быстрое разложение):\n";
    auto LU_dense_fast_C = LU_fast_C.toDense();
    for (const auto& row : LU_dense_fast_C) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "Результат умножения L * U (обычное разложение):\n";
    auto LU_dense_C = LU_C.toDense();
    for (const auto& row : LU_dense_C) {
        for (const auto& val : row) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(6) << val << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    C.toLIL();
    std::cout << "Матрицы C и L*U (быстрое разложение) "
              << (C == LU_fast_C ? "совпадают" : "не совпадают") << "\n";
    LU_fast_C.printDense(6);
    std::cout << "\n";

    std::cout << "Результат умножения L * U (обычное разложение):\n";
    LU_C.printDense(6);
    std::cout << "\n";


    std::cout << "Матрицы C и L*U (обычное разложение) " 
              << (C == LU_C ? "совпадают" : "не совпадают") << "\n";

    /*if (!instabilities_C.empty()) {
        std::cout << "\nОбнаружены места численной неустойчивости в быстром разложении:\n";
        for (const auto& inst : instabilities_C) {
            std::cout << "Шаг " << inst.k << ": "
                      << "Строка " << inst.i << ", "
                      << "Столбец " << inst.j << ", "
                      << "Значение = " << inst.value << ", "
                      << "Порог = " << inst.threshold << "\n";
        }
    }*/
}

int main() {
    // Установка русского языка в консоль
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    setlocale(LC_ALL, "Russian");

    testMatrix();

    return 0;
}