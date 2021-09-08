epsilon = 0.0000001


def run():
    matrix, frees = read_matrix("input.txt")
    u_matrix, l_matrix = decompose(matrix)

    print("Multiply: ")
    print_matrix(matrix_mul(l_matrix, u_matrix))

    solve_system(matrix, l_matrix, u_matrix, frees)


def read_matrix(filename):
    file = open(filename)

    matrix = []
    frees = []

    for line in file.readlines():
        row = list(map(lambda n: int(n), line.strip().split(" ")))
        matrix.append(row[:-1])
        frees.append(row[-1])

    print("Input matrix: ")
    print_matrix(matrix)

    print("Input frees: ")
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

    print("U-matrix: ")
    print_matrix(u_matrix)

    print("L-matrix: ")
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


def solve_system(matrix, l_matrix, u_matrix, frees):
    y_vector = perform_y(frees, l_matrix)


def perform_y(frees, l_matrix):
    y_vector = [i for i in frees]

    for i in range(len(y_vector)):
        for k in range(i):
            y_vector[i] -= l_matrix[i][k] * y_vector[k]

    print("Y vector: ")
    print(y_vector)

    return y_vector

if __name__ == '__main__':
    run()
