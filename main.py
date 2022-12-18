import math
import toml


def distance3d(vec1, vec2):
    return math.sqrt(((vec1[0] - vec2[0]) ** 2) + ((vec1[1] - vec2[1]) ** 2) + ((vec1[2] - vec2[2]) ** 2))


class Solver:
    def __init__(self, model_data):
        self.L = model_data["L"]
        self.L_k = model_data["L_k"]
        self.K = model_data["K"]
        self.K_l = model_data["K_l"]
        self.C = model_data["C"]
        self.a = model_data["a"]
        self.a_milled = [[i for _ in range(len(self.K))] for i in self.a]
        self.d = model_data["d"]

    def solve(self, accuracy=0.01, debug_info = False):
        b_start = [0 for i in range(len(self.K))]
        i = 1
        while True:
            print("\t===Итерация {}===".format(i))
            h_data = self.__h_calc()
            b_data = self.__B_calc(h_data)
            b_star_data = self.__B_star_calc(b_data)
            dist = distance3d(b_start, b_star_data)
            b_start = b_star_data
            print("h: {}".format(h_data))
            print("B_k(L): {}".format(b_data))
            print("a_k(L): {}".format(self.a_milled))
            print("B*: {}".format(b_star_data))
            print("\t===Конец итерации===\n")
            if dist <= accuracy:
                break
            self.a_milled = self.__a_milled_calc(b_data)
            i += 1
        return b_start

    def __h(self, n, l):
        if n < 0:
            return 0
        elif n == 0:
            return 1
        else:
            s = 0
            for i in self.K_l[l]:
                s += self.d[i] * self.a_milled[i][l] * self.__h(n - self.d[i], l)
            return (1.0 / n) * s

    def __h_calc(self):
        h_ = []
        for i in self.L:
            h_ln = []
            for j in range(self.C[i] + 1):
                h_ln.append(self.__h(j, i))
            h_.append(h_ln)
        return h_

    def __B(self, k, l, h_data):
        s1 = 0
        s2 = 0
        if self.C[l] - self.d[k] + 1 < 0:
            raise Exception("Bad arguments")
        for i in range(self.C[l] - self.d[k] + 1, self.C[l] + 1):
            s1 += h_data[l][i]
        for i in range(self.C[l] + 1):
            s2 += h_data[l][i]
        return s1 / s2

    def __B_calc(self, h_data):
        b_ = []
        for i in self.L:
            b_kl = []
            for j in self.K:
                b_kl.append(self.__B(i, j, h_data))
            b_.append(b_kl)
        return b_

    def __B_star_calc(self, b_data):
        b = []
        for i in self.K:
            p = 1
            for j in self.L_k[i]:
                p *= (1 - b_data[i][j])
            b.append(1 - p)
        return b

    def __a_milled_calc(self, b_data):
        a_milled_new = []
        for i in self.K:
            a_milled_new_small = []
            for j in self.L:
                p = self.a[i]
                for z in self.L_k[i]:
                    if z != j:
                        p *= (1 - b_data[i][z])
                a_milled_new_small.append(p)
            a_milled_new.append(a_milled_new_small)
        return a_milled_new


def main():
    tt = toml.load("model.toml")
    solve = Solver(tt["model"])
    blocks = solve.solve(tt["solver"]["accuracy"], tt["solver"]["debug_info"])

    print("Для данной системы: ")
    for i in range(len(blocks)):
        print("b_{} = {}".format(i, blocks[i]))


if __name__ == '__main__':
    main()
