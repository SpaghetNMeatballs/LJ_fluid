def aad(x_exp, x_calc):
    assert (len(x_exp) == len(x_calc))
    returnable = 0
    for i in range(len(x_exp)):
        returnable += abs((x_exp[i] - x_calc[i]))
        # / x_exp[i])
    returnable *= 100 / len(x_exp)
    return returnable
