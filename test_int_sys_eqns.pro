;test_int_sys_eqns.pro

;compile the various functions and procedures
.run odeint, odeint_derivs, int_sys_eqns

;name of procedure that calculates the derivatives dy/dx
pro_name = 'odeint_derivs'

;parameters used by pro_name are passed via this structure
a = 1.0
b = 2.0
c = 3.0
params = {a:a, b:b, c:c}

;initialize y(x=0)
y_init = fltarr(3)
y_init[0] = a
y_init[1] = 0.0
y_init[2] = c

;create x axis
N = 101
x_max = 10.0
x = x_max*findgen(N)/(N - 1)

;fractional error
error = 1.0e-4

;initial stepsize
dx_try = 1.0e-5

;minimum stepsize
dx_min = 1.0e-13

;Runga Kutta integration
int_sys_eqns, x, y_init, error, dx_try, dx_min, pro_name, y, params=params, quiet=0

;display results
window, xs=550, ys=850, retain=2
!p.multi = [0,1,3]

;expected solution y(x) = sin(x) + a
plot, x, y[0,*], charsize=2
oplot, x, y[0,*], thick=5
y_exp=sin(x) + a
oplot, x, y_exp, color=128, thick=2

;expected solution y(x) = b*x^2
plot, x, y[1,*], charsize=2
oplot, x, y[1,*], thick=5
y_exp=b*(x^2)
oplot, x, y_exp, color=128, thick=2

;expected solution y(x) = c*exp(-x)
plot, x, y[2,*], charsize=2
oplot, x, y[2,*], thick=5
y_exp=c*exp(-x)
oplot, x, y_exp, color=128, thick=2
!p.multi=0

