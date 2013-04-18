pro odeint_derivs, x, y, dydx, params=params

;extract parameters from params structure
a = params.a
b = params.b
c = params.c

;calculate derivatives dy/dx
dydx = fltarr(3)
dydx[0] = a*cos(x)
dydx[1] = 2.0*b*x
dydx[2] = -y[2]

return
end
