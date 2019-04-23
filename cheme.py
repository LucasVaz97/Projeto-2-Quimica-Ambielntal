from chempy import balance_stoichiometry
import numpy as np
import math as m


def mmc(num1, num2):
    mdc = num1
    y = num2

    resto = None
    while resto is not 0:
        resto = mdc % y
        mdc = y
        y = resto

    return (num1 * num2) / mdc


def parse(A):
    counter = 0
    if A[len(A)-1] == "+":
        return 1
    for i in A:
        if A[counter] == "+":
            return int(A[counter+1])
        counter = counter+1


class Metal:
    def __init__(self, name, oxi, pot):
        self.name = name
        self.oxi = oxi
        self.pot = pot


class Pilha:
    def __init__(self, metal1, metal2, Con1, Con2, temp, stack):
        self.stack = stack
        self.catodo = None
        self.anodo = None
        self.reativos = None
        self.produtos = None
        self.ddp = 0
        self.maxA = self.ddp/0.00001
        self.ConCatodo = None
        self.ConAnodo = None
        self.temp = temp+273.15
        self.massa = 1  # CONSERTAR
        IdxM1 = self.stack.index(metal1.name)
        IdxM2 = self.stack.index(metal2.name)
        if IdxM1 < IdxM2:
            self.catodo = metal2
            self.anodo = metal1
            self.ConCatodo = Con2
            self.ConAnodo = Con1
        else:
            self.catodo = metal1
            self.anodo = metal2
            self.ConCatodo = Con1
            self.ConAnodo = Con2

        print("Anodo:", self.anodo.name)
        print("Catodo:", self.catodo.name)
        self.reativos = {self.anodo.name, self.catodo.oxi}
        self.produtos = {self.anodo.oxi, self.catodo.name}
        print("Reagentes:", self.reativos)
        print("Produtos:", self.produtos)

        if(self.catodo.name != self.anodo.name):
            self.balance = self.Balance()

        self.DDP()

        self.carga = self.Carga()
        print("Carga: ", self.carga)
        self.capacidadeDeCarga = self.CalcCapacidadeDeCarga()
        self.densidadeDeCarga = self.CalcDensidadeDeCarga()
        self.densidadeDeEnergia = self.CalcDensidadeDeEnergia()

    def DDP(self):
        Const = (8.314*self.temp)/96500
        Eo = self.catodo.pot - self.anodo.pot
        print("DDP Teorica: ", Eo)
        n = mmc(parse(self.anodo.oxi), parse(self.catodo.oxi))
        if(self.anodo.name != self.catodo.name):
            anodo_exp = int(self.balance[1][self.anodo.oxi])
            catodo_exp = int(self.balance[1][self.catodo.name])
            log = np.log((self.ConAnodo**anodo_exp) /
                         (self.ConCatodo**catodo_exp))
        else:
            if(self.ConCatodo > self.ConAnodo):
                maiorC = self.ConCatodo
                menorC = self.ConAnodo
                maior_exp = self.catodo
                menor = self.anodo
            else:
                menorC = self.ConCatodo
                maiorC = self.ConAnodo
                maior = self.catodo
                menor = self.anodo

            log = np.log(menorC/maiorC)

        self.ddp = Eo-(Const/n)*log
        print("DDP REAL: ", self.ddp)

    def Carga(self):

        elec = mmc(parse(self.anodo.oxi), parse(self.catodo.oxi))
        const = 9.65e4
        columb = const * elec
        Ah = columb/3600

        return Ah

    def CalcCapacidadeDeCarga(self):
        return 0

    def CalcDensidadeDeCarga(self):
        return self.carga / self.massa

    def CalcDensidadeDeEnergia(self):
        return self.densidadeDeCarga * self.ddp

    def Custo(self):
        pass

    def Balance(self):
        reag, prod = balance_stoichiometry(self.reativos, self.produtos)
        return reag, prod


Zn = Metal('Zn', 'Zn+2', -0.76)
Cu = Metal('Cu', 'Cu+2', 0.34)
Pb = Metal('Pb', 'Pb+2', -0.13)
stack = ['Zn', 'Cr', 'Ni', 'Pb', 'Cu']
p1 = Pilha(Cu, Zn, 2, 2, 25, stack)
