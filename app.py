import random as rnd

epsilon = 0.0000001
is_debug = False
rounding = 5
file_out = open("out.txt", "w")


def run():
    file_in = open("input.txt")
    number_of_tests = int(file_in.readline())

    for _ in range(number_of_tests):
        matrix, frees = read_matrix(file_in)

        file_out.write("Исходная матрица: \n")
        print_matrix(matrix)
        file_out.write("Свободные члены: ")
        file_out.write(str(frees) + "\n")

        try:
            u_matrix, l_matrix = decompose(matrix)

            if is_debug:
                file_out.write("Невязка (исходная - полученная): \n")
                print_matrix(matrix_sub(matrix, matrix_mul(l_matrix, u_matrix)))

            solve_system(matrix, l_matrix, u_matrix, frees)
            determinant(u_matrix)
            matrix_reversed = matrix_reverse(l_matrix, u_matrix, matrix)
            norms(matrix, matrix_reversed)
        except Exception as e:
            pass

        file_out.write("\n****************************************\n\n")


def read_matrix(file):
    matrix = []
    frees = []

    number_of_rows = int(file.readline())

    if number_of_rows > 0:
        for i in range(number_of_rows):
            row = list(map(lambda n: float(n), file.readline().strip().split(" ")))
            matrix.append(row[:-1])
            frees.append(row[-1])

        if is_debug:
            file_out.write("Исходная матрица: \n")
            print_matrix(matrix)

            file_out.write("Свободные члены: \n")
            file_out.write(str(frees))
    else:
        # RANDOM MATRIX
        number_of_rows = int(file.readline())

        for i in range(number_of_rows):
            new_row = []
            for j in range(number_of_rows):
                new_row.append(round(rnd.randint(-100, 100) * rnd.random(), rounding))
            matrix.append(new_row)
            frees.append(round(rnd.randint(-100, 100) * rnd.random(), rounding))

    return matrix, frees


def print_matrix(matrix):
    file_out.write("\n".join(["\t".join(map(str, r)) for r in matrix]) + "\n")


def decompose(matrix):
    u_matrix, l_matrix = [[0 for __ in range(len(matrix))] for _ in range(len(matrix))], [
        [0 for __ in range(len(matrix))] for _ in range(len(matrix))]

    for i in range(len(matrix)):
        l_matrix[i][i] = 1

    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i <= j:
                u_matrix[i][j] = get_u(matrix, l_matrix, u_matrix, i, j)
                if i == j:
                    if abs(u_matrix[i][j]) < epsilon:
                        file_out.write("Один из главных миноров матрицы A равен 0, разложение невозможно\n")
                        raise Exception()
            else:
                l_matrix[i][j] = get_l(matrix, l_matrix, u_matrix, i, j)

    if is_debug:
        file_out.write("U-матрица: \n")
        print_matrix(u_matrix)

        file_out.write("L-матрица: \n")
        print_matrix(l_matrix)

    return u_matrix, l_matrix


def get_u(matrix, l_matrix, u_matrix, i, j):
    sum = 0
    for k in range(i):
        sum = round(sum + l_matrix[i][k] * u_matrix[k][j], rounding)

    sum = round(matrix[i][j] - sum, rounding)
    return sum


def get_l(matrix, l_matrix, u_matrix, i, j):
    sum = 0
    for k in range(j):
        sum = round(sum + l_matrix[i][k] * u_matrix[k][j], rounding)

    sum = round(matrix[i][j] - sum, rounding)
    sum = round(sum / u_matrix[j][j], rounding)

    return sum


def matrix_mul(m1, m2):
    res = []
    if len(m2) != len(m1[0]):
        file_out.write("Матрицы не могут быть перемножены\n")
        raise Exception()
    else:
        for z in range(0, len(m1)):
            t = []
            for j in range(len(m2[0])):
                s = 0
                for i in range(len(m1[0])):
                    s = round(s + m1[z][i] * m2[i][j], rounding)
                t.append(s)
            res.append(t)
    return res


def matrix_sub(m1, m2):
    m_res = []
    for i in range(len(m1)):
        m_res.append([])
        for j in range(len(m1[0])):
            m_res[i].append(round(m1[i][j] - m2[i][j], rounding))

    return m_res


def solve_system(matrix, l_matrix, u_matrix, frees):
    y_vector = perform_y(frees, l_matrix)
    x_vector = perform_x(y_vector, u_matrix)

    file_out.write("Решение СЛАУ: " + str(x_vector) + "\n")

    result_frees = matrix_mul(matrix, [[i] for i in x_vector])
    file_out.write("Невязки: " + str(matrix_sub(result_frees, [[i] for i in frees])) + "\n")


def perform_y(frees, l_matrix):
    y_vector = [0 for _ in frees]

    for i in range(len(y_vector)):
        sum = 0
        for k in range(i):
            sum = round(sum + l_matrix[i][k] * y_vector[k], rounding)
        y_vector[i] = round(frees[i] - sum, rounding)

    if is_debug:
        file_out.write("Y вектор: \n")
        file_out.write(str(y_vector) + "\n")

    return y_vector


def perform_x(y_vector, u_matrix):
    x_vector = [0 for _ in y_vector]

    for i in range(len(y_vector)):
        x_vector_index = len(y_vector) - 1 - i
        sum = 0
        for k in range(x_vector_index + 1, len(y_vector)):
            sum = round(sum + u_matrix[x_vector_index][k] * x_vector[k], rounding)
        x_vector[x_vector_index] = round((y_vector[x_vector_index] - sum) / u_matrix[x_vector_index][x_vector_index], rounding)

    return x_vector


def determinant(u_matrix):
    det = 1
    for i in range(len(u_matrix)):
        det = round(det * u_matrix[i][i], rounding)

    file_out.write("Определитель: " + str(det) + "\n")


def matrix_reverse(l_matrix, u_matrix, matrix):
    result = [[0 for _ in l_matrix[0]] for __ in l_matrix]

    for index in range(len(result)):
        i = len(result) - 1 - index

        perform_diagonal(i, u_matrix, result)

        for row in range(1, i + 1):
            perform_column(i - row, i, u_matrix, result)

        for col in range(1, i + 1):
            perform_row(i, i - col, l_matrix, result)

    file_out.write("Обратная матрица:\n")
    print_matrix(result)

    if is_debug:
        file_out.write("Невязки:\n")
        print_matrix(matrix_mul(matrix, result))

    return result


def perform_diagonal(j, u_matrix, result):
    result[j][j] = 1

    for k in range(j + 1, len(result)):
        result[j][j] = round(result[j][j] - u_matrix[j][k] * result[k][j], rounding)
    result[j][j] = round(result[j][j] / u_matrix[j][j], rounding)


def perform_column(i, j, u_matrix, result):
    result[i][j] = 0

    for k in range(i + 1, len(result)):
        result[i][j] = round(result[i][j] + u_matrix[i][k] * result[k][j], rounding)
    result[i][j] = round(-(result[i][j] / u_matrix[i][i]), rounding)


def perform_row(i, j, l_matrix, result):
    result[i][j] = 0

    for k in range(j + 1, len(result)):
        result[i][j] = round(result[i][j] - result[i][k] * l_matrix[k][j], rounding)


def norms(matrix, matrix_reversed):
    norm_matrix = max_norm(matrix)
    norm_matrix_reversed = max_norm(matrix_reversed)

    file_out.write("Нормы: \n")
    file_out.write("||А||: " + str(norm_matrix) + "\n")
    file_out.write("||А^-1||: " + str(norm_matrix_reversed) + "\n")
    file_out.write("||A|| * ||А^-1|| = " + str(norm_matrix * norm_matrix_reversed) + "\n")


def max_norm(matrix):
    max_vector = []
    for j in range(len(matrix[0])):
        col_sum = 0
        for i in range(len(matrix)):
            col_sum += matrix[i][j]
        max_vector.append(col_sum)

    return max(max_vector)


if __name__ == '__main__':
    run()
