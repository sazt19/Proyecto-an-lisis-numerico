from django.shortcuts import render
from django.http import HttpResponse
import sympy as sympy
from sympy import Symbol, Abs, sympify
from sympy.abc import x
import numpy as np
from scipy import linalg
# Create your views here.

def index(request):
    return HttpResponse("")

def biseccion(fx, Tol, Niter, a, b):
    output = {
        "columns": ["iter", "a", "xm", "b", "f(xm)", "E"],
        "iterations": Niter,
        "errors": list()
    }

    # Configuraciones iniciales
    datos = list()
    x = Symbol('x')
    i = 1
    error = 1.0000000
    Fun = sympify(fx)

    Fa = Fun.subs(x, a) # Funcion evaluada en a
    Fa = Fa.evalf()

    xm0 = 0.0
    Fxm = 0

    xm = (a + b)/2 # Punto intermedio

    Fxm = Fun.subs(x, xm) # Funcion evaluada en Xm
    Fxm = Fxm.evalf()

    try:
        datos.append([0, '{:^15.7f}'.format(a), '{:^15.7f}'.format(xm), '{:^15.7f}'.format(b), '{:^15.7E}'.format(Fxm)]) # Datos con formato dado
        while (error > Tol) and (i < Niter): # Se repite hasta que el intervalo sea lo pequeño que se desee
            if (Fa*Fxm < 0): # Se elecciona un intervalo inicial, donde el valor de la funcion cambie de signo en [a,b]
                b = xm
            else:
                a = xm # Cambia de signo en [m,b]

            xm0 = xm
            xm = (a+b)/2 # Se calcula el punto intermedio del intervalo - Divide el intervalo a la mitadd

            Fxm = Fun.subs(x, xm)
            Fxm = Fxm.evalf() # Se evalua el punto intermedio en la funcion

            error = Abs(xm-xm0) # Se calcula el error

            datos.append([i, '{:^15.7f}'.format(a), '{:^15.7f}'.format(xm), 
                            '{:^15.7f}'.format(b), '{:^15.7E}'.format(Fxm), '{:^15.7E}'.format(error)]) # Se van agregando las soluciones con el formato deseado

            i += 1
    except BaseException as e:
        if str(e) == "can't convert complex to float":
            output["errors"].append(
                "Error in data: found complex in calculations")
        else:
            output["errors"].append("Error in data: " + str(e))

        return output

    output["results"] = datos
    output["root"] = xm
    return output