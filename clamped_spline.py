# -*- coding: utf-8 -*-

import numpy as np

def get_clamped_spline_coef(x, y):
    if len(x) != len(y):
        return "The sizes of the input ranges are different"

    # Los parámetros x e y se transforman en numpy arrays.
    x = np.array(x)
    y = np.array(y)
    
    n = len(x) - 1;

    # Se definen las variables que se retornan
    a = np.zeros(n + 1)
    b = np.zeros(n + 1)
    c = np.zeros(n + 1)
    d = np.zeros(n + 1)

    # Se comienza
    fp1 = 0.0
    fpn = 0.0
    # Modificar vectores a, b, c
    # vector <double> h(n), l(n+1), u(n+1);
    h = np.zeros(n)
    l = np.zeros(n + 1)
    u = np.zeros(n + 1)
    
    # vector <double> alpha(n+1), z(n+1);
    alpha = np.zeros(n + 1)
    z = np.zeros(n + 1)

    fp1 = (y[1] - y[0]) / (x[1] - x[0])
    fpn = (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  
    # for (unsigned int i=0; i < n; i++)
        # h.at(i) = x[i + 1] - x[i];
    for i in range(0, n):
        h[i] = x[i + 1] - x[i]
        
    alpha[0] = ((y[1] - y[0]) * 3) / h[0] - fp1 * 3
    alpha[n] = fpn * 3 - ((y[n] - y[n - 1]) * 3) / h[n - 1]

    # for (unsigned int i=1; i < n ; i++)
    for i in range(1, n):
        alpha[i] = ((y[i + 1] - y[i])) * (3.0 / h[i]) - ((y[i] - y[i - 1])) * (3.0 / h[i - 1])                                  
        
    l[0] = 2 * h[0]
    u[0] = 0.5
    z[0] = alpha[0] / l[0]

    # for (unsigned int i=1; i < n ; i++)
    for i in range(1, n):
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1]
        u[i] = h[i] / l[i]
        z[i] = (alpha[i] - z[i - 1] * h[i - 1]) / l[i]

    l[n] = h[n - 1] * (2 - u[n - 1])
    z[n] = (alpha[n] - z[n - 1] * h[n - 1]) / l[n]
    c[n] = z[n]
        
    # for (j=1; j < (n+1) ; j++) 
    for j in range(1, n + 1):
        i = n - j
        c[i] = z[i] - c[i + 1] * u[i]
        b[i] = (y[i + 1] - y[i]) / h[i] - ((c[i + 1] + c[i] * 2) * h[i])/ 3
        d[i] = (c[i + 1] - c[i]) /(3 *h[i])
        a[i] = y[i]
    
    b[n] = d[n] = 0
    a[n] = y[n]
    
    return a, b, c, d
    
    
def get_clamped_spline_interpol(plazos, tasas):
    result = get_clamped_spline_coef(plazos, tasas)
    def cs(t):
        ind = np.where(plazos<=t)
        if len(ind[0]) == 0:
            index = 0
        else:
            index = np.max(ind)
        t -= plazos[index]
        valor = result[0][index] + result[1][index] * t + result[1][index] * t**2 + result[2][index] * t**3
        der = result[1][index] + 2 * result[1][index] * t + 3 * result[2][index] * t**2
        der2 = 2 * result[1][index] + 6 * result[2][index] * t
        return valor, der, der2
    return cs
