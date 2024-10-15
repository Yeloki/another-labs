import lab1, lab2, lab3, lab4

n, k = map(int, input("Введите номер ЛР и подномер задачи ЛР через пробел: ").split())

infile = input("Введите название файла с входными данными\n"
               "оставьте поле пустым, если\n"
               "хотите использовать дефолтное значение: ")

outfile = input("Введите название файла с входными данными\n"
               "оставьте поле пустым, если\n"
               "хотите использовать дефолтное значение (укажите 1 для использования консоли): ")

if n == 1:
    if k == 1:
        if not infile:
            infile = "lab11.in"
        if not outfile:
            outfile = "lab11.out"
        if outfile != '1':
            lab1.read_and_run11('input/' + infile, 'output/' + outfile)
        else:
            lab1.read_and_run11('input/' + infile, None)
    if k == 2:
        if not infile:
            infile = "lab12.in"
        if not outfile:
            outfile = "lab12.out"
        if outfile != '1':
            lab1.read_and_run12('input/' + infile, 'output/' + outfile)
        else:
            lab1.read_and_run12('input/' + infile, None)
    if k == 3:
        if not infile:
            infile = "lab13.in"
        if not outfile:
            outfile = "lab13.out"
        if outfile != '1':
            lab1.read_and_run13('input/' + infile, 'output/' + outfile)
        else:
            lab1.read_and_run13('input/' + infile, None)
    if k == 4:
        if not infile:
            infile = "lab14.in"
        if not outfile:
            outfile = "lab14.out"
        if outfile != '1':
            lab1.read_and_run14('input/' + infile, 'output/' + outfile)
        else:
            lab1.read_and_run14('input/' + infile, None)
    if k == 5:
        if not infile:
            infile = "lab15.in"
        if not outfile:
            outfile = "lab15.out"
        if outfile != '1':
            lab1.read_and_run15('input/' + infile, 'output/' + outfile)
        else:
            lab1.read_and_run15('input/' + infile, None)
