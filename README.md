# multiplication_matrix

Реализация алгоритма умножения матриц с блочным разделением данных с использованием MPI.

1. Реализовать алгоритм умножения матриц с блочным разделением данных
1. С++ с использованием MPI
1. Тест.

# Сборка

Сборка программы осуществляться средствами сборочной утилиты Cmake.

Для этого в каталоге с исходным кодом нужно запустить следующую команду:

```
cmake -DCMAKE_CXX_FLAGS:STRING="${CFLAGS}" ./
```

# Тестирование

В каталоге `test` находится программа `random_matrix.py` при запуске которой создаются файлы `A.txt`, 'B.txt' в которых находятся сгенерированные матрицы.

Также программа создаёт файл `C.txt`, содержимое которого можно сравнить с выводом программы основной программы.

# Выполнение

Для запуска программы необходимо ввести следующую команду, находясь в каталоге с исходным кодом:

```
mpiexec -np 4 ./multiplication_matrix ./test/A.txt ./test/B.txt
```
