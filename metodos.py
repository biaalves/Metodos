from sympy import *

y, t= symbols('y, t')



def Euler(y0, t0, h, qtde, func):
	vetor_x=[]
	vetor_y=[]
	vetor_x.append(t0)
	vetor_y.append(y0)
	y_aux=y0
	t_aux=t0	
	
	for i in range(0, qtde):
		y_aux= y_aux+h*func.subs([(y, y_aux), (t, t_aux)])
		t_aux=t_aux+h
		vetor_x.append(t_aux)
		vetor_y.append(y_aux)

	ref_saida = open("saida.txt", "a")		
	print("Metodo de Euler", file=ref_saida)
	print("Y(" + str(t0) +")=" +str(y0), file=ref_saida)
	print("h=" +str(h), file=ref_saida)
	
	for j in range(0, qtde+1):
		print(j, vetor_y[j], file=ref_saida)		
	ref_saida.write('\n')
	ref_saida.close()

def Euler_inverso(y0, t0, h, qtde, func):
	vetor_x=[]
	vetor_y=[]
	vetor_x.append(t0)
	vetor_y.append(y0)
	#for i in range(1, qtde):
	print("eae")	

def Euler_aprimorado(y0, t0, h, qtde, func):
	print("eae")

def Runge_Kutta(y0, t0, h, qtde, func):
	vetor_x=[]
	vetor_y=[]
	vetor_x.append(t0)
	vetor_y.append(y0)
	y_aux=y0
	t_aux=t0
	for i in range(0, qtde):
		k1= func.subs([(y, y_aux), (t, t_aux)])
		k2= func.subs([(y, y_aux+((h/2)*k1)), (t, t_aux+(h/2))])
		k3= func.subs([(y, y_aux+((h/2)*k2)), (t, t_aux+(h/2))])
		k4= func.subs([(y, y_aux+(h*k3)), (t, t_aux+h)])
	
		y_aux= y_aux + ((h/6)*(k1+(2*k2)+(2*k3)+k4))
		t_aux= t_aux+h
		vetor_x.append(t_aux)
		vetor_y.append(y_aux)

	ref_saida = open("saida.txt", "a")		
	print("Metodo de Runge Kutta", file=ref_saida)
	print("Y(" + str(t0) +")=" +str(y0), file=ref_saida)
	print("h=" +str(h), file=ref_saida)
	
	for j in range(0, qtde+1):
		print(j, vetor_y[j], file=ref_saida)	
	
	ref_saida.write('\n')
	ref_saida.close()
	
	
	
def main():	
	ref_arquivo = open("entrada.txt", "r")
	
	for linha in ref_arquivo:
		valores= linha.split()		
		if valores[0]== 'euler':
			Euler(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
		#elif valores[0] == 'euler_inverso':
		#	Euler_inverso(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
		#if valores[0] == 'euler_aprimorado':
		#	Euler_aprimorado(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
		elif valores[0] == 'runge_kutta':
			Runge_Kutta(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))

	#ref_arquivo.close()
if __name__ == '__main__': 
	main()

