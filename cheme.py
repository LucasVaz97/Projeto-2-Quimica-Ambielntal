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
    def __init__(self, name, oxi, pot, massaMolar):
        self.name = name
        self.oxi = oxi
        self.pot = pot
        self.massaMolar = massaMolar
        self.eletrons = parse(self.oxi)


class Pilha:
    def __init__(self, metal1, metal2, Con1, Con2, temp, massa1, massa2, stack):
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
        self.massaCloro = 35.453
        self.massaEletrodos = massa1 + massa2  # Eletrodos
        self.volumeSol = 0.001  # L
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

        self.massa = self.CalcMassaDaSolucao() + self.massaEletrodos
        self.reativos = {self.anodo.name, self.catodo.oxi}
        self.produtos = {self.anodo.oxi, self.catodo.name}

        if(self.catodo.name != self.anodo.name):
            self.balance = self.Balance()

        self.DDP()

        self.capacidadeDeCarga = self.CalcCapacidadeDeCarga()
        self.densidadeDeCarga = self.CalcDensidadeDeCarga()
        self.densidadeDeEnergia = self.CalcDensidadeDeEnergia()

        self.custo = self.Custo()

    def DDP(self):
        Const = (8.314*self.temp)/96500
        Eo = self.catodo.pot - self.anodo.pot
        n = mmc(self.anodo.eletrons, self.catodo.eletrons)
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

    def CalcCapacidadeDeCarga(self):

        elec = mmc(self.anodo.eletrons, self.catodo.eletrons)
        const = 9.65e4
        columb = const * elec
        Ah = columb/3600

        return 0

    def CalcDensidadeDeCarga(self):
        return self.capacidadeDeCarga / self.massa

    def CalcDensidadeDeEnergia(self):
        return self.densidadeDeCarga * self.ddp

    def Custo(self):
        return 0

    def Balance(self):
        reag, prod = balance_stoichiometry(self.reativos, self.produtos)
        return reag, prod

    def CalcMassaDaSolucao(self):
        m1 = self.volumeSol * self.ConAnodo * \
            (self.anodo.massaMolar + self.massaCloro * self.anodo.eletrons)
        m2 = self.volumeSol * self.ConCatodo * \
            (self.catodo.massaMolar + self.massaCloro * self.anodo.eletrons)
        return m1 + m2


def EscolhePilha(ddp, pot, temp):
    return 0


Li = Metal('Li', 'Li+', -3.04, 6.941)
K = Metal('K', 'K+', -2.92, 39.0983)
Ba = Metal('Ba', 'Ba+2', -2.9, 137.327)
Sr = Metal('Sr', 'Sr+2', -2.89, 87.62)
Ca = Metal('Ca', 'Ca+2', -2.87, 40.078)
Zn = Metal('Zn', 'Zn+2', -0.76, 65.38)
Cr = Metal('Cr', 'Cr+3', -0.74, 51.9961)
Fe = Metal('Fe', 'Fe+2', -0.44, 137.327)
Ni = Metal('Ni', 'Ni+2', -0.23, 58.6934)
Pb = Metal('Pb', 'Pb+2', -0.13, 207.2)
Cu = Metal('Cu', 'Cu+2', 0.34, 63.546)

stack = ['Li', 'K', 'Ba', 'Sr', 'Ca', 'Zn', 'Cr', 'Fe', 'Ni', 'Pb', 'Cu']
metais = {'Li': Li, 'K': K, 'Ba': Ba, 'Sr': Sr, 'Ca': Ca, 'Zn': Zn,
          'Cr': Cr, 'Fe': Fe, 'Ni': Ni, 'Pb': Pb, 'Cu': Cu}


def main():
    print("---------- Projeto Pilha ----------")
    print("")
    selection = int(input(
        "Deseja escolher uma pilha [0] ou saber as propriedades de uma pilha [1]? "))
    if selection == 1:
        metal1_name = input(
            "Escolha o primeiro metal ['Li', 'K', 'Ba', 'Sr', 'Ca', 'Zn', 'Cr', 'Fe', 'Ni', 'Pb', 'Cu']: ")
        metal1 = metais[metal1_name]
        con_metal1 = int(
            input("Qual a concentracao do primeiro metal em mol/L? "))
        massa1 = int(input("Qual a massa do primeiro eletrodo em Kg? "))
        metal2_name = input(
            "Escolha o segundo metal ['Li', 'K', 'Ba', 'Sr', 'Ca', 'Zn', 'Cr', 'Fe', 'Ni', 'Pb', 'Cu']: ")
        metal2 = metais[metal2_name]
        con_metal2 = int(
            input("Qual a concentracao do segundo metal em mol/L? "))
        massa2 = int(input("Qual a massa do segundo eletrodo em Kg? "))
        temperatura = int(
            input("Qual a temperatura de operacao em graus celsius? "))
        print("")

        print("---------- Detalhes da Pilha ----------")
        print("")
        pilha = Pilha(metal1, metal2, con_metal1,
                      con_metal2, temperatura, massa1, massa2, stack)
        print("Anodo:", pilha.anodo.name)
        print("Catodo:", pilha.catodo.name)
        print("Reagentes:", pilha.reativos)
        print("Produtos:", pilha.produtos)
        print("DDP: ", pilha.ddp)
        print("Corrente Maxima: ", pilha.maxA)
        print("Capacidade de Carga: ", pilha.capacidadeDeCarga)
        print("Densidade de Carga: ", pilha.densidadeDeCarga)
        print("Densidade de Energia: ", pilha.densidadeDeEnergia)
        print("Custo: ", pilha.custo)

    elif selection == 0:
        ddp = int(input("Qual a ddp desejada em V? "))
        pot = int(input("Qual a potencia desejada? "))
        temp = int(input("Quanto tempo em horas deve ficar ligado? "))
        pilha = EscolhePilha(ddp, pot, temp)


main()
