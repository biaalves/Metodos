from sympy import *

y, t = symbols('y, t')

arquivo_de_saida = None

def Euler(y0, t0, h, qtde, func):
	vetor_x=[]
	vetor_y=[]
	vetor_x.append(t0)
	vetor_y.append(y0)
	y_aux=y0
	t_aux=t0	
	
	for i in range(0, qtde):
		y_aux = y_aux + h * func.subs([(y, y_aux), (t, t_aux)])
		t_aux = t_aux+h
		vetor_x.append(t_aux)
		vetor_y.append(y_aux)
	
	print("Metodo de Euler", file=arquivo_de_saida)
	print("y(" + str(t0) + ")=" + str(y0), file=arquivo_de_saida)
	print("h =", h, file=arquivo_de_saida)
	
	for j in range(0, qtde+1):
		print(j, vetor_y[j], file=arquivo_de_saida)		
	arquivo_de_saida.write('\n')

def Euler_inverso(y0, t0, h, qtde, func):
	vetor_x=[]
	vetor_y=[]
	vetor_x.append(t0)
	vetor_y.append(y0)
	y_aux = y0
	t_aux = t0
	for i in range(0, qtde):
		# com metodo de previsao
		# y_prev = y_aux + func.subs([(y, y_aux), (t, t_aux)]) * h
		# t_aux = t_aux + h
		# y_aux = y_aux + func.subs([(y, y_prev), (t, t_aux)]) * h
		t_aux = t_aux + h
		yn1 = solve(y_aux + h * (func.subs(t, t_aux) ) - y, y)
		y_aux = yn1[0]
		vetor_y.append(y_aux)
		vetor_x.append(t_aux)

	print("Metodo de Euler Inverso", file=arquivo_de_saida)
	print("y(" + str(t0) + ") = " + str(y0), file=arquivo_de_saida)
	print("h = " + str(h), file=arquivo_de_saida)
	
	for j in range(0, qtde+1):
		print(j, vetor_y[j], file=arquivo_de_saida)		
	arquivo_de_saida.write('\n')

def Euler_aprimorado(y0, t0, h, qtde, func):
	vetor_x=[]
	vetor_y=[]
	vetor_x.append(t0)
	vetor_y.append(y0)
	y_aux = y0
	t_aux = t0
	for i in range(0, qtde):
		tn=t_aux
		t_aux = t_aux + h
		# com metodo de previsao
		y_prev = y_aux + func.subs([(y, y_aux), (t, tn)]) * h		
		y_aux = y_aux + (func.subs([(y, y_prev), (t, t_aux)])+func.subs([(y, y_aux), (t, tn)])) * h/2
		#yn1=solve(y_aux + h * (func.subs(t,t_aux) + func.subs([(t, tn), (y, y_aux)])) / 2 - y, y)
		#y_aux = yn1[0]
		vetor_y.append(y_aux)
		vetor_x.append(t_aux)

	print("Metodo de Euler Aprimorado", file=arquivo_de_saida)
	print("y(" + str(t0) + ") = " + str(y0), file=arquivo_de_saida)
	print("h = " + str(h), file=arquivo_de_saida)
	
	for j in range(0, qtde+1):
		print(j, vetor_y[j], file=arquivo_de_saida)		
	arquivo_de_saida.write('\n')

def Runge_Kutta(y0, t0, h, qtde, func):
	vetor_x=[]
	vetor_y=[]
	vetor_x.append(t0)
	vetor_y.append(y0)
	y_aux = y0
	t_aux = t0
	for i in range(0, qtde):
		k1= func.subs([(y, y_aux), (t, t_aux)])
		k2= func.subs([(y, y_aux+((h/2)*k1)), (t, t_aux+(h/2))])
		k3= func.subs([(y, y_aux+((h/2)*k2)), (t, t_aux+(h/2))])
		k4= func.subs([(y, y_aux+(h*k3)), (t, t_aux+h)])
	
		y_aux = y_aux + ((h/6)*(k1+(2*k2)+(2*k3)+k4))
		t_aux = t_aux + h
		vetor_x.append(t_aux)
		vetor_y.append(y_aux)
	
	print("Metodo de Runge Kutta", file=arquivo_de_saida)
	print("Y(" + str(t0) + ") = " + str(y0), file=arquivo_de_saida)
	print("h = " + str(h), file=arquivo_de_saida)
	
	for j in range(0, qtde+1):
		print(j, vetor_y[j], file=arquivo_de_saida)	
	
	arquivo_de_saida.write('\n')
	
def Adam_Bashforth(ys, t0, h, qtde, func, ordem):
	vetor_x=[t0]
	vetor_y=[]
	if ordem == 2:
		t_vetor=[t0-h, t0]		
		vetor_y.append(ys[1])
		for i in range(0, qtde):
			fn1 = 3/2 * (func.subs([(t, t_vetor[i+1]), (y, ys[i+1])]))
			fn= 1/2 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			yn2 = ys[i+1] + h * (fn1 - fn) 
			ys.append(yn2)
			tn=t_vetor[i+1] + h
			t_vetor.append(tn)
			vetor_y.append(yn2)
			vetor_x.append(tn)
	elif ordem == 3:
		t_vetor=[t0 - 2*h, t0-h, t0]		
		vetor_y.append(ys[2])
		for i in range(0, qtde):
			fn2 = 23/12 * (func.subs([(t, t_vetor[i+2]), (y, ys[i+2])]))
			fn1 = 4/3 * (func.subs([(t, t_vetor[i+1]), (y, ys[i+1])]))
			fn= 5/12 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			yn3 = ys[i+2] + h * (fn2 - fn1 + fn) 
			ys.append(yn3)
			tn=t_vetor[i+2] + h
			t_vetor.append(tn)
			vetor_y.append(yn3)
			vetor_x.append(tn)
		
	elif ordem == 4:
		t_vetor=[t0 - 3*h, t0 - 2*h, t0-h, t0]		
		vetor_y.append(ys[3])
		for i in range(0, qtde):
			fn3 = 55/24 * (func.subs([(t, t_vetor[i+3]), (y, ys[i+3])]))			
			fn2 = 59*/24 * (func.subs([(t, t_vetor[i+2]), (y, ys[i+2])]))
			fn1 = 37/24 * (func.subs([(t, t_vetor[i+1]), (y, ys[i+1])]))
			fn= 3/8 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			yn4 = ys[i+3] + h * (fn3 - fn2 + fn1 - fn) 
			ys.append(yn4)
			tn=t_vetor[i+3] + h
			t_vetor.append(tn)
			vetor_y.append(yn4)
			vetor_x.append(tn)

	elif ordem == 5:	
		t_vetor=[t0-4*h, t0-3*h, t0-2*h, t0-h, t0]		
		vetor_y.append(ys[4])
		
		for i in range(0, qtde):
			fn4 = 1901/720 * (func.subs([(t, t_vetor[i+4]), (y, ys[i+4])]))
			fn3 = 1387/360 * (func.subs([(t, t_vetor[i+3]), (y, ys[i+3])]))
			fn2 = 109/30 * (func.subs([(t, t_vetor[i+2]), (y, ys[i+2])]))
			fn1 = 637/360 * (func.subs([(t, t_vetor[i+1]), (y, ys[i+1])]))
			fn = 251/720 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			yn5 = ys[i+4] + h * (fn4 - fn3 + fn2 - fn1 + fn)
			ys.append(yn5)
			tn = t_vetor[i+4] + h
			t_vetor.append(tn)
			vetor_y.append(yn5)
			vetor_x.append(tn)			
	
	elif ordem == 6:
		t_vetor=[t0 - 5*h, t0-4*h, t0-3*h, t0-2*h, t0-h, t0]		
		vetor_y.append(ys[5])
		
		for i in range(0, qtde):
			fn5 = 4277/1440 * (funs.subs([(t, t_vetor[i+5]), (y, ys[i+5])]))
			fn4 = 2641/480 * (func.subs([(t, t_vetor[i+4]), (y, ys[i+4])]))
			fn3 = 4991/720 * (func.subs([(t, t_vetor[i+3]), (y, ys[i+3])]))
			fn2 = 3649/720 * (func.subs([(t, t_vetor[i+2]), (y, ys[i+2])]))
			fn1 = 959/480 * (func.subs([(t, t_vetor[i+1]), (y, ys[i+1])]))
			fn = 95/288 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			yn6 = ys[i+5] + h * (fn5 - fn4 + fn3 - fn2 + fn1 - fn)
			ys.append(yn6)
			tn = t_vetor[i+5] + h
			t_vetor.append(tn)
			vetor_y.append(yn6)
			vetor_x.append(tn)	

	elif ordem == 7:
		t_vetor = [t0 - 6*h, t0 - 5*h, t0-4*h, t0-3*h, t0-2*h, t0-h, t0]		
		vetor_y.append(ys[6])
		
		for i in range(0, qtde):
			fn6 = 198721/60480 * (funs.subs([(t, t_vetor[i+6]), (y, ys[i+6])]))
			fn5 = 18637/2520 * (funs.subs([(t, t_vetor[i+5]), (y, ys[i+5])]))
			fn4 = 235183/20160 * (func.subs([(t, t_vetor[i+4]), (y, ys[i+4])]))
			fn3 = 10754/945 * (func.subs([(t, t_vetor[i+3]), (y, ys[i+3])]))
			fn2 = 135713/20160 * (func.subs([(t, t_vetor[i+2]), (y, ys[i+2])]))
			fn1 = 5603/2520 * (func.subs([(t, t_vetor[i+1]), (y, ys[i+1])]))
			fn = 19087/60480 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			yn7 = ys[i+6] + h * (fn6 - fn5 + fn4 - fn3 + fn2 - fn1 + fn)
			ys.append(yn7)
			tn = t_vetor[i+6] + h
			t_vetor.append(tn)
			vetor_y.append(yn7)
			vetor_x.append(tn)	
	
	elif ordem == 8:
		t_vetor = [t0 - 7*h, t0 - 6*h, t0 - 5*h, t0-4*h, t0-3*h, t0-2*h, t0-h, t0]		
		vetor_y.append(ys[7])
		
		for i in range(0, qtde):
			fn7 = 16083/4480 * (funs.subs([(t, t_vetor[i+7]), (y, ys[i+7])]))
			fn6 = 1152169/120960 * (funs.subs([(t, t_vetor[i+6]), (y, ys[i+6])]))
			fn5 = 242653/13440 * (funs.subs([(t, t_vetor[i+5]), (y, ys[i+5])]))
			fn4 = 296053/13440 * (func.subs([(t, t_vetor[i+4]), (y, ys[i+4])]))
			fn3 = 2102243/120960 * (func.subs([(t, t_vetor[i+3]), (y, ys[i+3])]))
			fn2 = 115747/13440 * (func.subs([(t, t_vetor[i+2]), (y, ys[i+2])]))
			fn1 = 32863/13440 * (func.subs([(t, t_vetor[i+1]), (y, ys[i+1])]))
			fn = 5257/17280 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			yn8 = ys[i+7] + h * (fn7 - fn6 + fn5 - fn4 + fn3 - fn2 + fn1 - fn)
			ys.append(yn8)
			tn = t_vetor[i+7] + h
			t_vetor.append(tn)
			vetor_y.append(yn8)
			vetor_x.append(tn)	


	print("Metodo de Adam-Bashforth (ordem=", ordem, ")",  file=arquivo_de_saida)
	print("Y(" + str(t0) + ") = " + str(vetor_y[0]), file=arquivo_de_saida)
	print("h = " + str(h), file=arquivo_de_saida)
	
	for j in range(0, qtde+1):
		print(j, vetor_y[j], file=arquivo_de_saida)	
	
	arquivo_de_saida.write('\n')
	
def main():	
	global arquivo_de_saida
	arquivo_de_saida = open("saida.txt", "w")
	arquivo_de_entrada = open("entrada.txt", "r")
	
	for linha in arquivo_de_entrada:
		valores = linha.split()		
		if valores[0] == 'euler':
			Euler(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
		elif valores[0] == 'euler_inverso':
			Euler_inverso(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
		if valores[0] == 'euler_aprimorado':
			Euler_aprimorado(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
		elif valores[0] == 'runge_kutta':
			Runge_Kutta(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
		elif valores[0] == 'adam_bashforth':
			ordem = int(valores[-1])	
			ys= valores[1 : ordem+1]			
			ys=list(map(float, ys))
			Adam_Bashforth(ys, float(valores[ordem+1]), float(valores[ordem+2]), int(valores[ordem+3]), sympify(valores[ordem+4]), ordem)
	arquivo_de_saida.close()
	
if __name__ == '__main__': 
	main()