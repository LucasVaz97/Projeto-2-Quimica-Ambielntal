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
    def __init__(self, name, oxi, pot, massaMolar, preco):
        self.name = name
        self.oxi = oxi
        self.pot = pot
        self.massaMolar = massaMolar
        self.eletrons = parse(self.oxi)
        self.massa = 0
        self.preco = preco


Li = Metal('Li', 'Li+', -3.04, 6.941, 0.27)
K = Metal('K', 'K+', -2.92, 39.0983, 1)
Ba = Metal('Ba', 'Ba+2', -2.9, 137.327, 0.55)
Sr = Metal('Sr', 'Sr+2', -2.89, 87.62, 1)
Ca = Metal('Ca', 'Ca+2', -2.87, 40.078, 0.2)
Zn = Metal('Zn', 'Zn+2', -0.76, 65.38, 0.053)
Cr = Metal('Cr', 'Cr+3', -0.74, 51.9961, 0.32)
Fe = Metal('Fe', 'Fe+2', -0.44, 137.327, 0.072)
Ni = Metal('Ni', 'Ni+2', -0.23, 58.6934, 0.077)
Pb = Metal('Pb', 'Pb+2', -0.13, 207.2, 0.0245)
Cu = Metal('Cu', 'Cu+2', 0.34, 63.546, 0.0976)

stack = ['Li', 'K', 'Ba', 'Sr', 'Ca', 'Zn', 'Cr', 'Fe', 'Ni', 'Pb', 'Cu']
metais = {'Li': Li, 'K': K, 'Ba': Ba, 'Sr': Sr, 'Ca': Ca, 'Zn': Zn,
          'Cr': Cr, 'Fe': Fe, 'Ni': Ni, 'Pb': Pb, 'Cu': Cu}


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
        self.massa1 = massa1
        self.massa2 = massa2
        self.metal1 = metal1
        self.metal2 = metal2
        self.volumeSol = 0.001  # L
        IdxM1 = self.stack.index(metal1.name)
        IdxM2 = self.stack.index(metal2.name)

        if IdxM1 < IdxM2:
            self.catodo = metal2
            self.anodo = metal1
            self.ConCatodo = Con2
            self.ConAnodo = Con1
            self.catodo.massa = massa2
            self.anodo.massa = massa1
        else:
            self.catodo = metal1
            self.anodo = metal2
            self.ConCatodo = Con1
            self.ConAnodo = Con2
            self.catodo.massa = massa1
            self.anodo.massa = massa2

        self.massa = self.CalcMassaDaSolucao() + self.massa1 + self.massa2
        self.reativos = {self.anodo.name, self.catodo.oxi}
        self.produtos = {self.anodo.oxi, self.catodo.name}

        if(self.catodo.name != self.anodo.name):
            self.balance = self.Balance()

        self.DDP()
        self.CalcCapacidadeDeCarga()
        self.CalcDensidadeDeCarga()
        self.CalcDensidadeDeEnergia()
        self.Custo()

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
            else:
                menorC = self.ConCatodo
                maiorC = self.ConAnodo

            log = np.log(menorC/maiorC)

        self.ddp = Eo-(Const/n)*log

    def CalcCapacidadeDeCarga(self):
        MolAnodo = self.anodo.massa/self.anodo.massaMolar
        MolCatodo = self.catodo.massa/self.catodo.massaMolar
        # Qunatidade de mols de anodo
        if self.catodo.name != self.anodo.name:
            QTMA = int(self.balance[0][self.anodo.name]) * \
                MolCatodo/int(self.balance[0][self.catodo.oxi])
            if QTMA < MolAnodo:
                limitante = self.catodo
                limitanteName = self.catodo.oxi
            else:
                limitante = self.anodo
                limitanteName = self.anodo.name

            elec = mmc(self.anodo.eletrons, self.catodo.eletrons) / \
                int(self.balance[0][limitanteName]) * \
                limitante.massa/limitante.massaMolar

        else:
            if self.anodo.massa < self.catodo.massa:
                limitante = self.anodo
            else:
                limitante = self.catodo
            elec = mmc(self.anodo.eletrons, self.catodo.eletrons) * \
                limitante.massa/limitante.massaMolar

        const = 9.65e4
        columb = const * elec
        Ah = columb/3600

        self.capacidadeDeCarga = Ah

    def CalcDensidadeDeCarga(self):
        self.densidadeDeCarga = self.capacidadeDeCarga / self.massa

    def CalcDensidadeDeEnergia(self):
        self.densidadeDeEnergia = self.densidadeDeCarga * self.ddp

    def Custo(self):
        preco1 = self.massa1 * self.metal1.preco
        preco2 = self.massa2 * self.metal2.preco
        precoTot = preco1 + preco2
        self.custo = precoTot

    def Balance(self):
        reag, prod = balance_stoichiometry(self.reativos, self.produtos)
        return reag, prod

    def CalcMassaDaSolucao(self):
        m1 = self.volumeSol * self.ConAnodo * \
            (self.anodo.massaMolar + self.massaCloro * self.anodo.eletrons)
        m2 = self.volumeSol * self.ConCatodo * \
            (self.catodo.massaMolar + self.massaCloro * self.anodo.eletrons)
        return m1 + m2


def EscolhePilha(ddp, pot, tempo):

    con = 1
    temperatura = 25
    massa = 50

    amp = pot / ddp
    amph = amp * tempo

    menor_custo = 999999999999999999999999999999999
    menor_metal1 = ""
    menor_metal2 = ""

    for i in stack:

        metal1 = metais[i]

        for j in stack:

            metal2 = metais[j]

            pilha = Pilha(metal1, metal2, con, con,
                          temperatura, massa, massa, stack)

            soma = 0
            counter_ddp = 0
            while (ddp - soma >= 0.0000001):
                soma += pilha.ddp
                if soma == 0:
                    break
                counter_ddp += 1

            soma = 0
            counter_amph = 0
            while (amph - soma >= 0.0000001):
                soma += pilha.capacidadeDeCarga
                if soma == 0:
                    break
                counter_amph += 1

            nPilhas = counter_ddp * counter_amph

            custoTotal = nPilhas * pilha.custo

            if (custoTotal < menor_custo) and (custoTotal != 0):
                menor_custo = custoTotal
                menor_metal1 = metal1.name
                menor_metal2 = metal2.name

    return menor_custo, menor_metal1, menor_metal2


def main():
    print("---------- Projeto Pilha ----------")
    print("")
    selection = int(input(
        "Deseja escolher uma pilha [0] ou saber as propriedades de uma pilha [1]? "))
    if selection == 1:
        metal1_name = input(
            "Escolha o primeiro metal ['Li', 'K', 'Ba', 'Sr', 'Ca', 'Zn', 'Cr', 'Fe', 'Ni', 'Pb', 'Cu']: ")
        metal1 = metais[metal1_name]
        con_metal1 = float(
            input("Qual a concentracao do primeiro metal em mol/L? "))
        massa1 = float(input("Qual a massa do primeiro eletrodo em g? "))
        metal2_name = input(
            "Escolha o segundo metal ['Li', 'K', 'Ba', 'Sr', 'Ca', 'Zn', 'Cr', 'Fe', 'Ni', 'Pb', 'Cu']: ")
        metal2 = metais[metal2_name]
        con_metal2 = float(
            input("Qual a concentracao do segundo metal em mol/L? "))
        massa2 = float(input("Qual a massa do segundo eletrodo em g? "))
        temperatura = float(
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
        print("DDP em V: ", pilha.ddp)
        print("Capacidade de Carga em Ah: ", pilha.capacidadeDeCarga)
        print("Densidade de Carga: ", pilha.densidadeDeCarga)
        print("Densidade de Energia: ", pilha.densidadeDeEnergia)
        print("Custo em dolares: ", pilha.custo)

    elif selection == 0:
        ddp = float(input("Qual a ddp desejada em V? "))
        pot = float(input("Qual a potencia desejada em W? "))
        temp = float(
            input("Quanto tempo em horas deve ficar ligado em Horas? "))
        pilha = EscolhePilha(ddp, pot, temp)

        menor_metal1 = pilha[1]
        menor_metal2 = pilha[2]
        menor_custo = pilha[0]
        print("---------- Pilha Recomendada ----------")
        print("")
        print("Metal 1: ", menor_metal1)
        print("Metal 2: ", menor_metal2)
        print("Custo em dolares: ", menor_custo)


main()
