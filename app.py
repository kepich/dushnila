epsilon = 0.0000001
is_debug = False


def run():
    matrix, frees = read_matrix("input.txt")
    u_matrix, l_matrix = decompose(matrix)

    if is_debug:
        print("Невязка (исходная - полученная): ")
        print_matrix(matrix_sub(matrix, matrix_mul(l_matrix, u_matrix)))

    solve_system(matrix, l_matrix, u_matrix, frees)
    determinant(u_matrix)
    matrix_reverse(l_matrix, u_matrix)


def read_matrix(filename):
    file = open(filename)

    matrix = []
    frees = []

    for line in file.readlines():
        row = list(map(lambda n: float(n), line.strip().split(" ")))
        matrix.append(row[:-1])
        frees.append(row[-1])

    if is_debug:
        print("Исходная матрица: ")
        print_matrix(matrix)

        print("Свободные члены: ")
        print(frees)

    return matrix, frees


def print_matrix(matrix):
    print("\n".join(["\t".join(map(str, r)) for r in matrix]))


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
                        print("Один из главных миноров матрицы A равен 0, разложение невозможно")
                        exit(-1)
            else:
                l_matrix[i][j] = get_l(matrix, l_matrix, u_matrix, i, j)

    if is_debug:
        print("U-матрица: ")
        print_matrix(u_matrix)

        print("L-матрица: ")
        print_matrix(l_matrix)

    return u_matrix, l_matrix


def get_u(matrix, l_matrix, u_matrix, i, j):
    sum = 0
    for k in range(i):
        sum += l_matrix[i][k] * u_matrix[k][j]

    sum = matrix[i][j] - sum
    return sum


def get_l(matrix, l_matrix, u_matrix, i, j):
    sum = 0
    for k in range(j):
        sum += l_matrix[i][k] * u_matrix[k][j]

    sum = matrix[i][j] - sum
    sum /= u_matrix[j][j]

    return sum


def matrix_mul(m1, m2):
    res = []
    if len(m2) != len(m1[0]):
        print("Матрицы не могут быть перемножены")
    else:
        for z in range(0, len(m1)):
            t = []
            for j in range(len(m2[0])):
                s = 0
                for i in range(len(m1[0])):
                    s += m1[z][i] * m2[i][j]
                t.append(s)
            res.append(t)
    return res


def matrix_sub(m1, m2):
    m_res = []
    for i in range(len(m1)):
        m_res.append([])
        for j in range(len(m1[0])):
            m_res[i].append(m1[i][j] - m2[i][j])

    return m_res


def solve_system(matrix, l_matrix, u_matrix, frees):
    y_vector = perform_y(frees, l_matrix)
    x_vector = perform_x(y_vector, u_matrix)

    print("Решение СЛАУ: " + str(x_vector))

    result_frees = matrix_mul(matrix, [[i] for i in x_vector])
    print("Невязки: " + str(matrix_sub(result_frees, [[i] for i in frees])))


def perform_y(frees, l_matrix):
    y_vector = [i for i in frees]

    for i in range(len(y_vector)):
        for k in range(i):
            y_vector[i] -= l_matrix[i][k] * y_vector[k]

    if is_debug:
        print("Y вектор: ")
        print(y_vector)

    return y_vector


def perform_x(y_vector, u_matrix):
    x_vector = [i for i in y_vector]

    for i in range(len(y_vector)):
        x_vector_index = len(y_vector) - 1 - i
        for k in range(x_vector_index + 1, len(y_vector)):
            x_vector[x_vector_index] -= u_matrix[x_vector_index][k] * y_vector[k]
        x_vector[x_vector_index] /= u_matrix[x_vector_index][x_vector_index]

    return x_vector


def determinant(u_matrix):
    det = 1
    for i in range(len(u_matrix)):
        det *= u_matrix[i][i]

    print("Определитель: " + str(det))


def matrix_reverse(l_matrix, u_matrix):
    result = [[0 for _ in l_matrix[0]] for __ in l_matrix]

    for index in range(len(result)):
        i = len(result) - 1 - index

        perform_diagonal(i, u_matrix, result)

        for row in range(1, i):
            perform_column(i - row, i, u_matrix, result)

        for col in range(1, i):
            perform_row(i, i - col, l_matrix, result)




def perform_diagonal(j, u_matrix, result):
    result[j][j] = 1

    for k in range(j + 1, len(result)):
        result[j][j] -= u_matrix[j][k] * result[k][j]
    result[j][j] /= u_matrix[j][j]


def perform_column(i, j, u_matrix, result):
    result[i][j] = 0

    for k in range(i + 1, len(result)):
        result[i][j] += u_matrix[i][k] * result[k][j]
    result[i][j] /= - u_matrix[i][i]


def perform_row(i, j, l_matrix, result):
    result[i][j] = 0

    for k in range(j + 1, len(result)):
        result[i][j] -= result[i][k] * l_matrix[k][j]


if __name__ == '__main__':
    run()
