# Fast LU Decomposition

Реализация быстрого LU-разложения для разреженных матриц с использованием символьного анализа.

## Описание

Проект содержит реализацию двух алгоритмов LU-разложения:
1. Полное LU-разложение с символьным анализом
2. Быстрое LU-разложение с использованием сохраненной схемы операций

## Особенности

- Поддержка разреженных матриц
- Символьный анализ для оптимизации
- Обработка численно неустойчивых элементов
- Реализация в форматах LIL и CSR

## Тесты

Проект включает тесты для различных типов матриц:
1. Разреженная матрица с диагональным преобладанием
2. Разреженная матрица с несколькими ненулевыми элементами
3. Разреженная матрица с численно неустойчивыми элементами
4. Более заполненная разреженная матрица

## Требования

- C++ компилятор с поддержкой C++11 или выше
- Visual Studio (для Windows) 