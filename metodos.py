from sympy import *

y, t = symbols('y, t')

arquivo_de_saida = None

def Printar_Arquivo(strr, y0, t0, h, qtde, vetor_y):
	print("Metodo de", strr, file=arquivo_de_saida)
	print("y(" + str(t0) + ") = " + str(y0), file=arquivo_de_saida)
	print("h = " + str(h), file=arquivo_de_saida)
	for j in range(0, qtde+1):
		print(j, vetor_y[j], file=arquivo_de_saida)		
	arquivo_de_saida.write('\n')


def Euler(y0, t0, h, qtde, func, cond):
	#cond eh a variavel que indica se a funcao foi chamada de adams(==1)
	vetor_x=[t0]
	vetor_y=[y0]	
	y_aux=y0
	t_aux=t0	
	
	for i in range(0, qtde):
		y_aux = y_aux + h * func.subs([(y, y_aux), (t, t_aux)])
		t_aux = t_aux+h
		
		vetor_x.append(t_aux)
		vetor_y.append(y_aux)
	if cond == 0: #nao veio de adams, printa no arquivo
		Printar_Arquivo('Euler', y0, t0, h, qtde, vetor_y)
	elif cond == 1: #retorna a lista de ys
		return vetor_y
def Euler_inverso(y0, t0, h, qtde, func, cond):
	vetor_x=[]
	vetor_y=[]
	vetor_x.append(t0)
	vetor_y.append(y0)
	y_aux = y0
	t_aux = t0
	for i in range(0, qtde):
		#com metodo de previsao
		#y_prev = y_aux + func.subs([(y, y_aux), (t, t_aux)]) * h
		#t_aux = t_aux + h
		#y_aux = y_aux + func.subs([(y, y_prev), (t, t_aux)]) * h
		t_aux = t_aux + h
		yn1 = solve(y_aux + h * (func.subs(t, t_aux) ) - y, y)
		y_aux = yn1[0]
		vetor_y.append(y_aux)
		vetor_x.append(t_aux)
	if cond == 1:
		return vetor_y
	elif cond == 0:
		Printar_Arquivo('Euler Inverso', y0, t0, h, qtde, vetor_y)

def Euler_Aprimorado(y0, t0, h, qtde, func, cond):
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
	if cond == 1:
		return vetor_y
	elif cond == 0:
		Printar_Arquivo('Euler Aprimorado', y0, t0, h, qtde, vetor_y)

def Runge_Kutta(y0, t0, h, qtde, func, cond):
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

	if cond == 1:
		return vetor_y
	elif cond == 0:
		Printar_Arquivo('Runge Kutta', y0, t0, h, qtde, vetor_y)
	
def Adam_Bashforth(ys, t0, h, qtde, func, ordem, cond):
	
	if ordem == 2:
		t_vetor=[t0, t0 + h]		
		
		for i in range(ordem-1, qtde):
			fn1 = 3/2 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			fn= 1/2 * (func.subs([(t, t_vetor[i-1]), (y, ys[i-1])]))
			yn2 = ys[i] + h * (fn1 - fn) 
			ys.append(yn2)
			tn=t_vetor[i] + h
			t_vetor.append(tn)
			
	elif ordem == 3:
		t_vetor=[t0, t0 + h, t0 + 2*h]		
		
		for i in range(ordem-1, qtde):
			fn2 = 23/12 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			fn1 = 4/3 * (func.subs([(t, t_vetor[i-1]), (y, ys[i-1])]))
			fn= 5/12 * (func.subs([(t, t_vetor[i-2]), (y, ys[i-2])]))
			yn3 = ys[i] + h * (fn2 - fn1 + fn) 
			ys.append(yn3)
			tn=t_vetor[i] + h
			t_vetor.append(tn)
					
	elif ordem == 4:
		t_vetor=[t0, t0 + h, t0 + 2*h, t0 + 3*h]		
		
		for i in range(ordem-1, qtde):
			fn3 = 55/24 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))			
			fn2 = 59/24 * (func.subs([(t, t_vetor[i-1]), (y, ys[i-1])]))
			fn1 = 37/24 * (func.subs([(t, t_vetor[i-2]), (y, ys[i-2])]))
			fn= 3/8 * (func.subs([(t, t_vetor[i-3]), (y, ys[i-3])]))
			yn4 = ys[i] + h * (fn3 - fn2 + fn1 - fn) 
			ys.append(yn4)
			tn=t_vetor[i] + h
			t_vetor.append(tn)			

	elif ordem == 5:	
		t_vetor=[t0, t0 + h, t0 + 2*h, t0 + 3*h, t0 + 4*h]			
		for i in range(ordem-1, qtde):
			fn4 = 1901/720 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			fn3 = 1387/360 * (func.subs([(t, t_vetor[i-1]), (y, ys[i-1])]))
			fn2 = 109/30 * (func.subs([(t, t_vetor[i-2]), (y, ys[i-2])]))
			fn1 = 637/360 * (func.subs([(t, t_vetor[i-3]), (y, ys[i-3])]))
			fn = 251/720 * (func.subs([(t, t_vetor[i-4]), (y, ys[i-4])]))
			yn5 = ys[i] + h * (fn4 - fn3 + fn2 - fn1 + fn)
			ys.append(yn5)
			tn = t_vetor[i] + h
			t_vetor.append(tn)
				
	elif ordem == 6:
		
		t_vetor=[t0, t0 + h, t0 + 2*h, t0 + 3*h, t0 + 4*h, t0 + 5*h]		
		
		for i in range(ordem-1, qtde):
			fn5 = (4277/1440) * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			fn4 = (2641/480) * (func.subs([(t, t_vetor[i-1]), (y, ys[i-1])]))
			fn3 = (4991/720) * (func.subs([(t, t_vetor[i-2]), (y, ys[i-2])]))
			fn2 = (3649/720) * (func.subs([(t, t_vetor[i-3]), (y, ys[i-3])]))
			fn1 = (959/480) * (func.subs([(t, t_vetor[i-4]), (y, ys[i-4])]))
			fn = (95/288) * (func.subs([(t, t_vetor[i-5]), (y, ys[i-5])]))
			yn6 = ys[i] + (h * (fn5 - fn4 + fn3 - fn2 + fn1 - fn))
			ys.append(yn6)
			tn = t_vetor[i] + h
			
			t_vetor.append(tn)
			

	elif ordem == 7:
		t_vetor = [t0, t0 + h, t0 + 2*h, t0 + 3*h, t0 + 4*h, t0 + 5*h, t0 + 6*h]		
			
		for i in range(ordem-1, qtde):
			fn6 = 198721/60480 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			fn5 = 18637/2520 * (func.subs([(t, t_vetor[i-1]), (y, ys[i-1])]))
			fn4 = 235183/20160 * (func.subs([(t, t_vetor[i-2]), (y, ys[i-2])]))
			fn3 = 10754/945 * (func.subs([(t, t_vetor[i-3]), (y, ys[i-3])]))
			fn2 = 135713/20160 * (func.subs([(t, t_vetor[i-4]), (y, ys[i-4])]))
			fn1 = 5603/2520 * (func.subs([(t, t_vetor[i-5]), (y, ys[i-5])]))
			fn = 19087/60480 * (func.subs([(t, t_vetor[i-6]), (y, ys[i-6])]))
			yn7 = ys[i+6] + h * (fn6 - fn5 + fn4 - fn3 + fn2 - fn1 + fn)
			ys.append(yn7)
			tn = t_vetor[i] + h
			t_vetor.append(tn)
				
	elif ordem == 8:
		t_vetor = [t0, t0 + h, t0 + 2*h, t0 + 3*h, t0 + 4*h, t0 + 5*h, t0 + 6*h, t0 + 7*h]		
				
		for i in range(0, qtde):
			fn7 = 16083/4480 * (func.subs([(t, t_vetor[i]), (y, ys[i])]))
			fn6 = 1152169/120960 * (func.subs([(t, t_vetor[i-1]), (y, ys[i-1])]))
			fn5 = 242653/13440 * (func.subs([(t, t_vetor[i-2]), (y, ys[i-2])]))
			fn4 = 296053/13440 * (func.subs([(t, t_vetor[i-3]), (y, ys[i-3])]))
			fn3 = 2102243/120960 * (func.subs([(t, t_vetor[i-4]), (y, ys[i-4])]))
			fn2 = 115747/13440 * (func.subs([(t, t_vetor[i-5]), (y, ys[i-5])]))
			fn1 = 32863/13440 * (func.subs([(t, t_vetor[i-6]), (y, ys[i-6])]))
			fn = 5257/17280 * (func.subs([(t, t_vetor[i-7]), (y, ys[i-7])]))
			yn8 = ys[i] + h * (fn7 - fn6 + fn5 - fn4 + fn3 - fn2 + fn1 - fn)
			ys.append(yn8)
			tn = t_vetor[i] + h
			t_vetor.append(tn)				

	if cond == 1:		
		return ys

	elif cond == 0:
		strr = 'Adam Bashforth (ordem = '+str(ordem)+')'
		Printar_Arquivo(strr, ys[0], t0, h, qtde, ys)
		

def Adam_Bashforth_by_method(y0, t0, h, qtde, func, ordem, strr):
	
	if strr == 'Euler':
		ys = Euler(y0, t0, h, ordem-1, func, 1)
	elif strr == 'Euler Inverso':
		ys = Euler_Inverso(y0, t0, h, ordem-1, func, 1)
	elif strr == 'Euler Aprimorado':
		ys = Euler_Aprimorado(y0, t0, h, ordem-1, func, 1)
	elif strr == 'Runge Kutta':
		ys = Runge_Kutta(y0, t0, h, ordem-1, func, 1)

	vetor_y = Adam_Bashforth(ys, t0, h, qtde, func, ordem, 1)
	strr_2= 'Adam Bashforth por ' + strr + '(ordem=' +str(ordem)+')'
	Printar_Arquivo(strr_2, vetor_y[0], t0, h, qtde, vetor_y)

def Adam_Multon(ys, t0, h, qtde, func, ordem, cond):
	if ordem == 2:
		t_vetor = [t0]
		for i in range(ordem-2, qtde):
			t_aux= t_vetor[i] + h #tn+1
			yn2 = solve( ys[i] + (h/2) * (func.subs(t, t_aux) + func.subs([(t, t_vetor[i]), (y, ys[i])]))-y, y)
			y_aux= yn2[0]
			ys.append(y_aux)
			t_vetor.append(t_aux)

	elif ordem == 3:
		t_vetor = [t0, t0+h]
		for i in range(ordem-2, qtde):
			t_aux=t_vetor[i] + h
			fn1 = 2/3 * func.subs([(t, t_vetor[i]), (y, ys[i])])
			fn = 1/12 * func.subs([(t, t_vetor[i-1]), (y, ys[i-1])])
			yn3 = solve( ys[i] + h * ( 5/12 * func.subs(t, t_aux) + fn1 - fn)-y, y)
			y_aux= yn3[0]
			ys.append(y_aux)
			t_vetor.append(t_aux)

	elif ordem == 4:
		t_vetor = [t0, t0+h, t0+2*h]
		for i in range(ordem-2, qtde):
			t_aux=t_vetor[i] + h
			fn2 = 19/24 * func.subs([(t, t_vetor[i]), (y, ys[i])])
			fn1 = 5/24 * func.subs([(t, t_vetor[i-1]), (y, ys[i-1])])
			fn = 1/24 * func.subs([(t, t_vetor[i-2]), (y, ys[i-2])])
			yn4 = solve( ys[i] + h * ( 3/8 * func.subs(t, t_aux) + fn2 - fn1 + fn)-y, y)
			y_aux= yn4[0]
			ys. append(y_aux)
			t_vetor.append(t_aux)
	elif ordem == 5:
		t_vetor = [t0, t0+h, t0+2*h, t0+3*h]
		for i in range(ordem-2, qtde):
			t_aux=t_vetor[i] + h
			fn3 = 323/360 * func.subs([(t, t_vetor[i]), (y, ys[i])])
			fn2 = 11/30 * func.subs([(t, t_vetor[i-1]), (y, ys[i-1])])
			fn1 = 53/360 * func.subs([(t, t_vetor[i-2]), (y, ys[i-2])])
			fn = 19/720 * func.subs([(t, t_vetor[i-3]), (y, ys[i-3])])
			yn5 = solve( ys[i] + h * ( 251/720 * func.subs(t, t_aux) + fn3 - fn2 + fn1 - fn)-y, y)
			y_aux= yn5[0]
			ys. append(y_aux)
			t_vetor.append(t_aux)
	elif ordem == 6:
		t_vetor = [t0, t0+h, t0+2*h, t0+3*h, t0+4*h]
		for i in range(ordem-2, qtde):
			t_aux=t_vetor[i] + h			
			fn4 = 1427/1440 * func.subs([(t, t_vetor[i]), (y, ys[i])])
			fn3 = 133/240 * func.subs([(t, t_vetor[i-1]), (y, ys[i-1])])
			fn2 = 241/720 * func.subs([(t, t_vetor[i-2]), (y, ys[i-2])])
			fn1 = 173/1440 * func.subs([(t, t_vetor[i-3]), (y, ys[i-3])])
			fn = 3/160 * func.subs([(t, t_vetor[i-4]), (y, ys[i-4])])
			yn6 = solve( ys[i] + h * ( 95/288 * func.subs(t, t_aux) + fn4 - fn3 + fn2 - fn1 + fn)-y, y)
			y_aux= yn6[0]
			ys. append(y_aux)
			t_vetor.append(t_aux)
	elif ordem == 7:
		t_vetor = [t0, t0+h, t0+2*h, t0+3*h, t0+4*h, t0+5*h]
		for i in range(ordem-2, qtde):
			t_aux=t_vetor[i] + h
			fn5 = 2713/2520 * func.subs([(t, t_vetor[i]), (y, ys[i])])
			fn4 = 15487/20160 * func.subs([(t, t_vetor[i-1]), (y, ys[i-1])])
			fn3 = 586/945 * func.subs([(t, t_vetor[i-2]), (y, ys[i-2])])
			fn2 = 6737/20160 * func.subs([(t, t_vetor[i-3]), (y, ys[i-3])])
			fn1 = 263/2520 * func.subs([(t, t_vetor[i-4]), (y, ys[i-4])])
			fn = 863/60480 * func.subs([(t, t_vetor[i-5]), (y, ys[i-5])])
			yn7 = solve( ys[i] + h * ( 19087/60480 * func.subs(t, t_aux) + fn5 - fn4 + fn3 - fn2 + fn1 - fn)-y, y)
			y_aux= yn7[0]
			ys. append(y_aux)
			t_vetor.append(t_aux)
	elif ordem == 8:
		t_vetor = [t0, t0+h, t0+2*h, t0+3*h, t0+4*h, t0+5*h, t0+6*h]
		for i in range(ordem-2, qtde):
			t_aux=t_vetor[i] + h
			fn6 = 139849/120960 * func.subs([(t, t_vetor[i]), (y, ys[i])])
			fn5 = 4511/4480 * func.subs([(t, t_vetor[i-1]), (y, ys[i-1])])
			fn4 = 123133/120960 * func.subs([(t, t_vetor[i-2]), (y, ys[i-2])])
			fn3 = 88547/120960 * func.subs([(t, t_vetor[i-3]), (y, ys[i-3])])
			fn2 = 1537/4480 * func.subs([(t, t_vetor[i-4]), (y, ys[i-4])])
			fn1 = 11351/120960 * func.subs([(t, t_vetor[i-5]), (y, ys[i-5])])
			fn = 275/24192 * func.subs([(t, t_vetor[i-6]), (y, ys[i-6])])
			yn8 = solve( ys[i] + h * ( 5257/17280 * func.subs(t, t_aux) + fn6 - fn5 + fn4 - fn3 + fn2 - fn1 + fn)-y, y)
			y_aux= yn8[0]
			ys. append(y_aux)
			t_vetor.append(t_aux)
	
	if cond == 0:
		strr = 'Adam Multon (ordem = '+str(ordem)+')'
		Printar_Arquivo(strr, ys[0], t0, h, qtde, ys)
	elif cond == 1:
		return ys

def Adam_Multon_by_method(y0, t0, h, qtde, func, ordem, strr):
	
	if strr == 'Euler':
		ys = Euler(y0, t0, h, ordem-2, func, 1)
	elif strr == 'Euler Inverso':
		ys = Euler_Inverso(y0, t0, h, ordem-2, func, 1)
	elif strr == 'Euler Aprimorado':
		ys = Euler_Aprimorado(y0, t0, h, ordem-2, func, 1)
	elif strr == 'Runge Kutta':
		ys = Runge_Kutta(y0, t0, h, ordem-2, func, 1)

	vetor_y = Adam_Multon(ys, t0, h, qtde, func, ordem, 1)
	strr_2= 'Adam Multon por ' + strr + '(ordem=' +str(ordem)+')'
	Printar_Arquivo(strr_2, vetor_y[0], t0, h, qtde, vetor_y)

def Formula_Inversa(ys, t0, h, qtde, func, ordem, cond):
	if ordem == 2:
		pass
	elif ordem == 3:
		pass
	elif ordem == 4:
		pass
	elif ordem == 5:
		pass
	elif ordem == 6:
		pass

#-------------------------MAIN---------------------------------
def main():	
	global arquivo_de_saida
	arquivo_de_saida = open("saida.txt", "w")
	arquivo_de_entrada = open("entrada.txt", "r")
	
	for linha in arquivo_de_entrada:
		valores = linha.split()		
		if valores[0] == 'euler':

			Euler(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), 0)
		elif valores[0] == 'euler_inverso':

			Euler_inverso(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), 0)
		elif valores[0] == 'euler_aprimorado':

			Euler_Aprimorado(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), 0)

		elif valores[0] == 'runge_kutta':

			Runge_Kutta(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), 0)

		elif valores[0] == 'adam_bashforth':

			ordem = int(valores[-1])	
			ys= valores[1 : ordem+1]			
			ys=list(map(float, ys))
			Adam_Bashforth(ys, float(valores[ordem+1]), float(valores[ordem+2]), int(valores[ordem+3]), sympify(valores[ordem+4]), ordem, 0)

		elif valores[0] == 'adam_bashforth_by_euler':		
	
			Adam_Bashforth_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), 'Euler')

		elif valores[0] == 'adam_bashforth_by_euler_inverso':		
	
			Adam_Bashforth_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), 'Euler Inverso')
		elif valores[0] == 'adam_bashforth_by_euler_aprimorado':		
	
			Adam_Bashforth_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), 'Euler Aprimorado')
		elif valores[0] == 'adam_bashforth_by_runge_kutta':		
	
			Adam_Bashforth_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), 'Runge Kutta')
	
		elif valores[0] == 'adam_multon':
			ordem = int(valores[-1])	
			ys= valores[1 : ordem]			
			ys=list(map(float, ys))			
			Adam_Multon(ys, float(valores[ordem]), float(valores[ordem+1]), int(valores[ordem+2]), sympify(valores[ordem+3]), ordem, 0)

		elif valores[0] == 'adam_multon_by_euler':		
	
			Adam_Multon_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), 'Euler')

		elif valores[0] == 'adam_multon_by_euler_inverso':		
	
			Adam_Multon_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), 'Euler Inverso')

		elif valores[0] == 'adam_multon_by_euler1_aprimorado':		
	
			Adam_Multon_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), 'Euler Aprimorado')

		elif valores[0] == 'adam_multon_by_runge_kutta':		
	
			Adam_Multon_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), 'Runge Kutta')
	
		elif valores[0] == 'formula_inversa':
			ordem = int(valores[-1])	
			ys= valores[1 : ordem]			
			ys=list(map(float, ys))			
			Formula_Inversa(ys, float(valores[ordem]), float(valores[ordem+1]), int(valores[ordem+2]), sympify(valores[ordem+3]), ordem, 0)
	arquivo_de_saida.close()
	
if __name__ == '__main__': 
	main()
