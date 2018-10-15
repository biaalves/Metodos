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
	if ordem == 2:
		
	elif ordem == 3:
	elif ordem == 4:
	elif ordem == 5:
	elif ordem == 6:

	
def main():	
	global arquivo_de_saida
	arquivo_de_saida = open("saida.txt", "w")
	ref_arquivo = open("entrada.txt", "r")
	
	for linha in ref_arquivo:
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
			ordem = valores[-1]
			ys= valores[1:ordem+1]
			ys=list(map(float, ys))
			Adam_Bashforth(ys, float(valores[ordem+1]), float(valores[oredem+2]), int(valores[ordem+3]), sympify(valores[ordem+4]), ordem)
	arquivo_de_saida.close()
	
if __name__ == '__main__': 
	main()

