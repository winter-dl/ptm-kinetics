def dotSeparation(txt1, txt2, *txts, numOfDots=60):
    numOfDots = numOfDots - len(txt1)
    print(txt1, end=" ")
    print(f"{numOfDots*'.'}", end=" ")
    if txts:
        print(txt2, end=" ")
        for txt in txts:
            print(txt, end=" ")
        print(f"\n", end="")
    else:
        print(txt2)

def x_for_maxabs_y(xs, func):
    ys = []
    for x in xs:
        y = func(x)
        y = abs(y)
        ys.append(y)
    x_y = dict(zip(ys, xs))
    fit = x_y[max(ys)]
    return fit

def x_for_minabs_y(xs, func):
    ys = []
    for x in xs:
        y = func(x)
        y = abs(y)
        ys.append(y)
    x_y = dict(zip(ys, xs))
    fit = x_y[min(ys)]
    return fit