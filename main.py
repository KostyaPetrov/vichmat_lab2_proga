import numpy as np
import matplotlib.pyplot as plt
# C://Users//Kostya//Desktop//ITMO//ITMO_4sem//VichMath//file.txt

def drow_graph(func, min_x, max_x, min_y, max_y, step):
    x = np.linspace(min_x, max_x, 10000)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.plot(x, func(x), "g", linewidth=2.0)

    ax.set(xlim=(min_x, max_x), xticks=np.arange(min_x, max_x, step),
           ylim=(min_y, max_y), yticks=np.arange(min_y, max_y, step))

    plt.show()


def input_interval(func):
    a, b = 0, 0
    correct = False
    while correct == False:
        try:
            a = float(input("Введите нижнюю границу интервала: "))
            b = float(input("Введите верхнюю границу интервала: "))
            if a > b:
                a, b = b, a
            elif a == b:
                raise ArithmeticError
            elif func(a) * func(b) >= 0:
                raise AttributeError
            correct = True
        except ValueError:
            print("Границы интервала должны быть числами")
        except ArithmeticError:
            print("Границы интервала не могут быть равны")
        except AttributeError:
            print("Интервал не содержит одного корня")
    return a, b


def find_proisv(func):
    return lambda x: (func(x + 0.000000001) - func(x - 0.000000001)) / (2 * 0.000000001)


def method_secushchih(func, a, b, acc):
    x_prev_prev = 0
    x_prev = 0
    # Поиск начального приближения
    if (func(a) * find_proisv(find_proisv(func))(a) > 0):
        print("В качестве начального приближения выбрано число: " + str(a))
        x_prev_prev = a
        x_prev = a + 0.001
    else:
        print("В качестве начального приближения выбрано число: " + str(b))
        x_prev_prev = b
        x_prev = b - 0.001

    x_current = x_prev * 1000 + 10
    x_prev = b
    iterations = 0

    # Поиск корней, отлавливаем условие остановки
    while (abs((func(x_prev))) > acc):
        x_current = x_prev - ((x_prev - x_prev_prev) / (func(x_prev) - func(x_prev_prev))) * func(x_prev)
        x_prev_prev = x_prev
        x_prev = x_current
        iterations += 1
    return x_current, iterations


def method_hord(func, a, b, acc):
    k = func(a) * find_proisv(find_proisv(func))(a)

    if (k > 0):
        c = a
        d = b
    else:
        c = b
        d = a
    interation = 0
    acc_cur = 1
    while acc_cur >= acc:
        x = c - func(c) * (d - c) / (func(d) - func(c))
        acc_cur = abs(x - c)
        c = x
        interation += 1
    return x, interation


def method_iterations(func, a, b, acc):
    fp_a = find_proisv(func)(a)
    fp_b = find_proisv(func)(b)

    print("Производная в точке A: " + str(fp_a))
    print("Производная в точке B: " + str(fp_b))

    lyambd_a = -(1 / fp_a)
    lyambd_b = - (1 / fp_b)

    print("Лямбда А = " + str(lyambd_a))
    print("Лямбда B = " + str(lyambd_b))

    if (fp_a > fp_b):
        lyambd = lyambd_a
    else:
        lyambd = lyambd_b

    fi = lambda x: x + lyambd * func(x)

    fi_s = find_proisv(fi)

    fi_s_a = fi_s(a)
    fi_s_b = fi_s(b)

    print("Производная фи в А: " + str(fi_s_a))
    print("Производная фи в B: " + str(fi_s_b))

    if (abs(fi_s(a)) > 1 or abs(fi_s(b)) > 1):
        print("Не удовлетворяет достаточному условию сходимости")
    else:
        print("Удовлетворяет достаточному условию сходимости")

    x_current = a
    x_prev = a * 1000 + 10

    iterations = 0

    # Поиск корней
    while (abs(x_prev - x_current) > acc) or (abs(func(x_current)) > acc):

        x_prev = x_current
        x_current = x_prev + lyambd * func(x_prev)
        iterations += 1
        if (iterations > 1000):
            print("Алгоритм расходится")
            exit()

    return x_current, iterations


# Для системы
# Частная производная по х
def find_dx(function, x, y, h=0.00000001):
    return (function(x + h, y) - function(x - h, y)) / (2 * h)


# Частная производная по у
def find_dy(function, x, y, h=0.00000001):
    return (function(x, y + h) - function(x, y - h)) / (2 * h)


# Поиск определителя матрицы Якоби
def opred(function1, function2, x, y):
    return find_dx(function1, x, y) * find_dy(function2, x, y) - find_dx(function2, x, y) * find_dy(function1, x, y)


def system_Nuthon(num_first_function, num_second_function, start_x, start_y, acc):
    x_current = start_x
    y_current = start_y
    y_prev = y_current * 1000 + 10
    x_prev = x_current * 1000 + 10
    iterations = 0

    while max(abs(x_current - x_prev), abs(y_current - y_prev)) > acc:
        x_prev = x_current
        y_prev = y_current
        J = opred(num_first_function, num_second_function, x_prev, y_prev)
        A = num_first_function(x_prev, y_prev) / J
        B = num_second_function(x_prev, y_prev) / J
        x_current = x_prev - A * find_dy(num_second_function, x_prev, y_prev) + B * find_dy(num_first_function, x_prev, y_prev)
        y_current = y_prev + A * find_dx(num_second_function, x_prev, y_prev) - B * find_dy(num_first_function, x_prev, y_prev)
        iterations += 1
        if (iterations == 100):
            print("Расходится")
            exit()
    return x_current, y_current, abs(x_current - x_prev), abs(y_current - y_prev), iterations


command = 0
while command != 1 and command != 2:
    try:
        command = int(input("Введите 1, чтобы выбрать одно уравнение. Введите 2, чтобы выбрать систему уравнений: "))
    except ValueError:
        print("Вы должны ввести 1 или 2")

if (command == 1):
    # Выбор уравнения
    print("Выберите уравнение из списка:")
    print("1 : 2,3*x^3 + 5,75*x^2 − 7,41*x − 10,6")
    print("2 : 2*x^2 - 5 * x - 4")
    print("3 : sin(x) - cos(x) + 0.2*x")

    num_func = -1
    while num_func != 1 and num_func != 2 and num_func != 3:

        try:
            num_func = int(input("Введите номер выбранной функции: "))
        except ValueError:
            print("Вы должны выбрать номер функции и ввести число")

    if (num_func == 1):
        func = lambda x: 2.3*x**3 + 5.75*x**2 - 7.41*x - 10.6
    elif (num_func == 2):
        func = lambda x: 2*x ** 2 - 5 * x - 4
    else:
        func = lambda x: np.sin(x) - np.cos(x) + 0.2*x

    # Выбор места чтения
    place_read_command = 0
    while place_read_command != 1 and place_read_command != 2:
        try:
            place_read_command = int(
                input("Введите 1, чтобы ввести данные с клавиатуры. Введите 2, чтобы прочитать данные из файла: "))
        except ValueError:
            print("Вы должны ввести 1 или 2")


    if(place_read_command==1):
        # Выбор способа выбора данных
        interval_command = 0
        while interval_command != 1 and interval_command != 2:
            try:
                interval_command = int(
                    input("Введите 1, чтобы ввести интервал. Введите 2, чтобы дать программе самой найти интервалы: "))
            except ValueError:
                print("Вы должны ввести 1 или 2")

        if (interval_command == 1):
            a, b = input_interval(func)
        else:
            corrrect = False

            for i in np.arange(-100, 100, 0.5):
                if (func(i) * func(i - 0.5) < 0):
                    print("Найден интервал: [" + str(i - 0.5) + ", " + str(i) + "]")
                    a = i - 0.5
                    b = i
                    corrrect = True
                    break

            if (corrrect == False):
                print("Интервал не найден")
                a, b = input_interval(func)

        #   Ввод точности
        accuracy = 0
        correct = False
        while correct == False:
            try:
                accuracy = float(input("Введите точность: "))
                if accuracy > 0:
                    correct = True
                else:
                    print("Число должно быть положительным")
            except ValueError:
                print("Введите число")
    else:
        correct=False
        while correct==False:
            file=str(input("Введите путь к файлу:"))
            with open(file, mode="rt") as fd:
                try:
                    a_str=fd.readline()
                    a_str = a_str.replace(",", ".")
                    a=float(a_str)
                    b_str=fd.readline()
                    b_str=b_str.replace(",", ".")
                    b=float(b_str)
                    accuracy=float(fd.readline())
                    if accuracy > 0:
                        correct = True
                    else:
                        print("Число должно быть положительным")
                        raise ArithmeticError
                    if a > b:
                        a, b = b, a
                    elif a == b:
                        raise ArithmeticError
                    elif func(a) * func(b) >= 0:
                        raise AttributeError
                    correct = True
                except ValueError:
                    print("Неверные данные в файле")
                except ArithmeticError:
                    print("Повторите попытку")

    # Выбор метода
    method = 0
    while method != 1 and method != 2 and method != 3:
        try:
            method = int(input(
                "Введите 1, чтобы выбрать метод хорд. Введите 2 чтобы выбрать метод секущих. Введите 3 чтобы выбрать метод простой итерации: "))
        except ValueError:
            print("Введите 1, 2 или 3")

    if (method == 1):
        answer, count_iterations = method_hord(func, a, b, accuracy)
    elif (method == 2):
        answer, count_iterations = method_secushchih(func, a, b, accuracy)
    else:
        answer, count_iterations = method_iterations(func, a, b, accuracy)

    print("Корень: " + str(answer) + " найден за " + str(count_iterations) + " итераций, f(x) = " + str(func(answer)))

    drow_graph(func, a - 0.5, b + 0.5, func(a) - 0.5, func(b) + 0.5, 0.2)
else:
    # Выбор уравнений для системы
    correct = False
    while correct == False:
        # Выбор уравнения
        print("Выберите уравнение из списка:")
        print("1 : x^3 + 2.28x^2 - 1.934x - 3.907")
        print("2 : x^2 - 3x - 2")
        print("3 : sin(x) - cos(x) + 0.2x")

        try:
            num_first_function = int(input("Введите номер первого уравнения: "))
            num_second_function = int(input("Введите номер второго уравнения: "))
            if (num_first_function == 1 or num_first_function == 2 or num_first_function == 3):
                if (num_second_function == 1 or num_second_function == 2 or num_second_function == 3):
                    if (num_first_function != num_second_function):
                        correct = True
                        break
            print("Числа должны быть разными и от 1 до 3")

        except ValueError:
            print("Должен быть числом")

    #   Получение начального приближения

    correct = False
    while correct == False:
        try:
            start_x = float(input("Введите начальное приближение: "))
            start_y = float()
            correct = True
        except ValueError:
            print("Начальное приближение должно быть числами, введенными через пробел.")

    #   Ввод точности
    accuracy = 0
    correct = False
    while correct == False:
        try:
            accuracy = float(input("Введите точность: "))
            if accuracy > 0:
                correct = True
            else:
                print("Число должно быть положительным")
        except ValueError:
            print("Введите число")

    x, y, acc_x, acc_y, count_iteration = system_Nuthon(num_first_function,num_second_function, start_x, start_y,accuracy)

    print("x = " + str(x) + ", y = " + str(y) + ", найден за " + str(count_iteration) + " итераций, вектор погрешностей: [" + str(acc_x) + ", " + str(acc_y) + "]")