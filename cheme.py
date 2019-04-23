from chempy import balance_stoichiometry
import numpy as np
import math as m
def mmc(num1, num2):
    mdc = num1
    y = num2

    resto = None
    while resto is not 0:
        resto = mdc % y
        mdc  = y
        y  = resto

    return (num1 * num2) / mdc


def parse(A):
    counter=0
    if A[len(A)-1]=="+":
        return 1
    for i in A:
        if A[counter]=="+":
            return int(A[counter+1])
        counter= counter+1

class Metal:
    def __init__(self,name,oxi,pot):
        self.name =name
        self.oxi=oxi
        self.pot=pot


class Pilha:
    def __init__(self,metal1,metal2,Con1,Con2,temp):
        self.stack=['Zn','Cr','Ni','Pb','Cu']
        self.catodo=None
        self.anodo=None
        self.reativos=0
        self.produtos=0
        self.ddp=0
        self.maxA=self.ddp/0.00001
        self.ConCatodo=None
        self.ConAnodo=None
        self.temp=temp+273.15
        IdxM1=self.stack.index(metal1.name)
        IdxM2=self.stack.index(metal2.name)
        if IdxM1<IdxM2:
            self.catodo=metal2
            self.anodo=metal1
            self.ConCatodo=Con2
            self.ConAnodo=Con1
        else:
            self.catodo=metal1
            self.anodo=metal2
            self.ConCatodo=Con1
            self.ConAnodo=Con2


    def DDP(self):
        maiorP=0
        menorP=0
        maiorC=0
        menorC=0
        print("Anodo:",self.anodo.name)
        print("Catodo:",self.catodo.name)
        Const=(8.314*self.temp)/96500
        self.reativos={self.anodo.name,self.catodo.oxi}
        self.produtos={self.anodo.oxi,self.catodo.name}
        print("Reagentes:",self.reativos)
        print("Produtos:",self.produtos)
        if(self.anodo.pot>self.catodo.pot):
            maiorP=self.anodo.pot
            menorP=self.catodo.pot
        else:
            menorP=self.anodo.pot
            maiorP=self.catodo.pot
        Eo=maiorP-menorP
        n=mmc(parse(self.anodo.oxi),parse(self.catodo.oxi))
        if(self.anodo.name != self.catodo.name):
            self.ddp=Eo-(Const/n)*np.log(self.ConAnodo/self.ConCatodo)
            print(Eo, self.ddp)
        else:
            if(self.ConCatodo>self.ConAnodo):
                maiorC=self.ConCatodo
                menorC=self.ConAnodo
            else:
                menorC=self.ConCatodo
                maiorC=self.ConAnodo
            self.ddp=Eo-(Const/n)*np.log(menorC/maiorC)
            print(Eo, self.ddp)
            

    def Carga(self):

        elec=parse(self.anodo.oxi)
        const = 9.65e4
        columb = const * elec
        Ah= columb/3600

        return Ah
        

            
    def Balance(self):
        reag,prod=balance_stoichiometry(self.reativos,self.produtos)
        return reag , prod
    





Zn=Metal('Zn','Zn+2',-0.76)
Cu=Metal('Cu','Cu+2',0.34)
Pb=Metal('Pb','Pb+2',-0.13)
p1=Pilha(Pb,Cu,2,2,25)
p1.DDP()
print(p1.Carga())




    


