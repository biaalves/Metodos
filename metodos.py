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

class MetodosPassoSimples:
  def euler(y0, t0, h, qtde, func):
      vetor_x=[t0]
      vetor_y=[y0]	
      y_aux=y0
      t_aux=t0	

      for i in range(0, qtde):
          y_aux = y_aux + h * func.subs([(y, y_aux), (t, t_aux)])
          t_aux = t_aux+h

          vetor_x.append(t_aux)
          vetor_y.append(y_aux)

      return vetor_y

  def euler_inverso(y0, t0, h, qtde, func):
      vetor_x=[t0]
      vetor_y=[y0]	
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

      return vetor_y

  def euler_aprimorado(y0, t0, h, qtde, func):
      vetor_x=[t0]
      vetor_y=[y0]	
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

      return vetor_y
	
  def runge_kutta(y0, t0, h, qtde, func):
      vetor_x=[t0]
      vetor_y=[y0]	
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


      return vetor_y
    
coeficientes_adam_bashforth = [
	[],
	[1],
	[3/2,          -1/2],
	[23/12,        -4/3,            5/12],
	[55/24,        -59/24,          37/24,        -3/8],
	[1901/720,     -1287/360,       109/30,       -637/360,      251/720],
	[4277/1440,    -2641/480,       4991/720,     -3649/720,     959/480,        -95/288],
	[198721/60480, -18637/2520,     235183/20160, -10754/945,    135713/20160,   -5603/2520,    19087/60480],
	[16083/4480,   -1152169/120960, 242653/13440, -296053/13440, 2102243/120960, -115747/13440, 32863/13440, -5257/17280]
]

def Adam_Bashforth(ys, t0, h, qtde, func, ordem):
	cur_coef = coeficientes_adam_bashforth[ordem] # pega os coeficients da ordem atual
	t_vetor = [t0] # cria o t_vetor inicial
	for i in range(1, ordem):
		t_vetor.append(t_vetor[i - 1] + h)
		
	for i in range(ordem-1, qtde):
		sum_fn = 0
		for j in range(ordem): # calcula a soma dos fn * coeficiente
			sum_fn += cur_coef[j] * func.subs([(t, t_vetor[i - j]), (y, ys[i - j])])
		yn = ys[i] + h * sum_fn # calcula novo y
		tn = t_vetor[i] + h # calcula novo t
		ys.append(yn)
		t_vetor.append(tn)
	return ys		

def Adam_Bashforth_by_method(y0, t0, h, qtde, func, ordem, strr):
	ys = getattr(MetodosPassoSimples, strr)(y0, t0, h, ordem-1, func)
	vetor_y = Adam_Bashforth(ys, t0, h, qtde, func, ordem)
	strr_2= 'Adam Bashforth por ' + strr + '(ordem=' +str(ordem)+')'
	Printar_Arquivo(strr_2, vetor_y[0], t0, h, qtde, vetor_y)

coeficientes_adam_moulton = [
	[],
	[1],
	[1/2,         1/2],
	[5/12,        2/3,           -1/12],
	[3/8,         19/24,         -5/24,        1/24],
	[251/720,     323/360,       -11/30,       53/360,        -19/720],
	[95/288,      1427/1440,     -133/240,     241/720,       -173/1440,     3/160],
	[19087/60480, 2713/2520,     -15487/20160, 586/945,       -6737/20160,   263/2520,  -863/60480],
	[5257/17280,  139849/120960, -4511/4480,   123133/120960, -88547/120960, 1537/4480, -11351/120960, 275/24192]
]

def Adam_Multon(ys, t0, h, qtde, func, ordem):
	cur_coef = coeficientes_adam_moulton[ordem] # pega os coeficientes da ordem atual
	t_vetor = [t0] # cria o t_vetor inicial
	for i in range(1, ordem - 1):
		t_vetor.append(t_vetor[i - 1] + h)
		
	for i in range(ordem-2, qtde):
		sum_fn = 0
		for j in range(ordem-1): # calcula a soma dos fn * coeficiente
			sum_fn += cur_coef[j + 1] * func.subs([(t, t_vetor[i - j]), (y, ys[i - j])])
		t_aux = t_vetor[i] + h
		yn = solve(ys[i] + h * (sum_fn + cur_coef[0] * func.subs(t, t_aux)) - y, y)# calcula novo y
		y_aux = yn[0]
		ys.append(y_aux)
		t_vetor.append(t_aux)
	return ys
  
def Adam_Multon_by_method(y0, t0, h, qtde, func, ordem, strr):
	ys = getattr(MetodosPassoSimples, strr)(y0, t0, h, ordem-2, func)
	vetor_y = Adam_Multon(ys, t0, h, qtde, func, ordem)
	strr_2= 'Adam Multon por ' + strr + '(ordem=' +str(ordem)+')'
	Printar_Arquivo(strr_2, vetor_y[0], t0, h, qtde, vetor_y)

coeficientes_formula_inversa = [
	[],
	[],
	[1,           1],
	[2/3,         4/3,            -1/3],
	[6/11,        18/11,          -9/11,        2/11],
	[12/25,       48/25,          -36/25,       16/25,      -3/25],
	[60/137,      300/137,        -300/137,     200/137,    -75/137,   12/137],
	[60/147,      360/147,        -450/147,     400/147,    -225/147,  72/147,    -10/147],
]
    

def Formula_Inversa(ys, t0, h, qtde, func, ordem):
	cur_coef = coeficientes_formula_inversa[ordem] # pega os coeficientes da ordem atual
	t_vetor = [t0] # cria o t_vetor inicial
	for i in range(1, ordem - 1):
		t_vetor.append(t_vetor[i - 1] + h)
		
	for i in range(ordem-2, qtde):
		sum_fn = 0
		for j in range(ordem-1): # calcula a soma dos yn * coeficiente
			sum_fn += cur_coef[j + 1] * ys[i-j]
		t_aux = t_vetor[i] + h
		yn = solve(sum_fn + h * ( cur_coef[0] * func.subs(t, t_aux)) - y, y)# calcula novo y
		y_aux = yn[0]
		ys.append(y_aux)
		t_vetor.append(t_aux)
	return ys

def Formula_Inversa_by_method(y0, t0, h, qtde, func, ordem, strr):
	ys = getattr(MetodosPassoSimples, strr)(y0, t0, h, ordem-2, func)
	vetor_y = Formula_Inversa(ys, t0, h, qtde, func, ordem)
	strr_2= 'Formula Inversa de Diferenciacao por ' + strr + '(ordem=' +str(ordem)+')'
	Printar_Arquivo(strr_2, vetor_y[0], t0, h, qtde, vetor_y)
#-------------------------MAIN---------------------------------
def main():	
	global arquivo_de_saida
	arquivo_de_saida = open("saida.txt", "w")
	arquivo_de_entrada = open("entrada.txt", "r")
	
	for linha in arquivo_de_entrada:
		valores = linha.split()		
		if valores[0] == 'euler':

			vetor_y = MetodosPassoSimples.euler(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
			Printar_Arquivo('Euler', float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), vetor_y)
		elif valores[0] == 'euler_inverso':

			vetor_y = MetodosPassoSimples.euler_inverso(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
			Printar_Arquivo('Euler Inverso', float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), vetor_y)
		elif valores[0] == 'euler_aprimorado':
			
			vetor_y = MetodosPassoSimples.euler_aprimorado(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
			Printar_Arquivo('Euler Aprimorado', float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), vetor_y)

		elif valores[0] == 'runge_kutta':

			vetor_y = MetodosPassoSimples.runge_kutta(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]))
			Printar_Arquivo('Runge Kutta', float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), vetor_y)

		elif valores[0] == 'adam_bashforth':

			ordem = int(valores[-1])	
			ys= valores[1 : ordem+1]			
			ys=list(map(float, ys))
			vetor_y = Adam_Bashforth(ys, float(valores[ordem+1]), float(valores[ordem+2]), int(valores[ordem+3]), sympify(valores[ordem+4]), ordem)
			strr = 'Adam Bashforth (ordem = '+str(ordem)+')'
			Printar_Arquivo(strr, vetor_y[0], float(valores[ordem+1]), float(valores[ordem+2]), int(valores[ordem+3]), vetor_y)

		elif valores[0].startswith('adam_bashforth_by_'):		
			metodo_auxiliar = valores[0][18:]
			Adam_Bashforth_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), metodo_auxiliar)		
	
		elif valores[0] == 'adam_multon':
			ordem = int(valores[-1])	
			ys = valores[1 : ordem]			
			ys = list(map(float, ys))			
			vetor_y = Adam_Multon(ys, float(valores[ordem]), float(valores[ordem+1]), int(valores[ordem+2]), sympify(valores[ordem+3]), ordem)
			strr = 'Adam Moulton (ordem = '+str(ordem)+')'
			Printar_Arquivo(strr, vetor_y[0], float(valores[ordem]), float(valores[ordem+1]), int(valores[ordem+2]), vetor_y)

		elif valores[0].startswith('adam_multon_by_'):		
			metodo_auxiliar = valores[0][15:]
			Adam_Multon_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), metodo_auxiliar)		
	
		elif valores[0] == 'formula_inversa':
			ordem = int(valores[-1])	
			ys= valores[1 : ordem]			
			ys=list(map(float, ys))			
			vetor_y = Formula_Inversa(ys, float(valores[ordem]), float(valores[ordem+1]), int(valores[ordem+2]), sympify(valores[ordem+3]), ordem)
			strr = 'Formula Inversa de Diferenciacao (ordem = '+str(ordem)+')'
			Printar_Arquivo(strr, vetor_y[0], float(valores[ordem]), float(valores[ordem+1]), int(valores[ordem+2]), vetor_y)
            
		elif valores[0].startswith('formula_inversa_by_'):		
			metodo_auxiliar = valores[0][19:]
			vetor_y = Formula_Inversa_by_method(float(valores[1]), float(valores[2]), float(valores[3]), int(valores[4]), sympify(valores[5]), int(valores[6]), metodo_auxiliar)		
		
	arquivo_de_saida.close()
	
if __name__ == '__main__': 
	main()
