;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Perform a Runga Kutta integration of a system of equations
;y=[y0(x), y1(x), y2(x),...], and sample the output at various times x. The inputs
;are x (which is a 1D vector where yi(x) is sampled, the boundary conditions
;y0 (which is a 1D vector that are the initial yi(x[0]) values), the allowed
;fractional error, the suggested initial timestep dx_try, the minimum timestep dx_min,
;pro_name=string that provides the name of the proceedure that calculates the 
;derivative dyi/dx (a 1D vector), and the 2D output array y[j,k], where j
;indicates the quantity yi which is evaluated at x[k]. Additional parameters
;are passes on to pro_name via params, which is often a structure. Set optional
;quiet=1 to avoid printing the code's %progress.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro int_sys_eqns, x, y0, error, dx_try, dx_min, pro_name, y, params=params, quiet=quiet

;create storage for y[j, x] array
Nx = n_elements(x)
Ny = n_elements(y0)
ytype = size(y0, type=1)
y = make_array(Ny, Nx, type=ytype)
fraction = 0

;store initial conditions
y[0, 0] = y0

;loop over all x and advance y[j] from x[j-1] to x[j]
for j = 1, Nx - 1 do begin $
    ystart = y[*, j - 1] & $
    x1 = x[j - 1] & $
    x2 = x[j] & $
    odeint, ystart, x1, x2, error, dx_try, dx_min, pro_name, params=params & $
    y[0, j] = ystart & $
    if (keyword_set(quiet) eq 0) then begin $
        fraction_old = fraction & $
        fraction = round(10.0*(x2 - x[0])/(x[Nx - 1] - x[0])) & $
        if (fraction ne fraction_old) then $
            print, format = '($, A0)', strtrim(string(10*fraction),2)+'% ' & $
    endif & $
endfor
if (keyword_set(quiet) eq 0) then print, '100%'

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
