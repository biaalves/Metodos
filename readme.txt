1. Antes de rodar o programa:
	- É necessário que esteja em ambiente Linux.

	- Versão do python: python 3.

	- Você deve ter no diretório do projeto os arquivos requirements.txt e RUNME
		* Navegue no terminal até o diretório do projeto e digite o comando- source RUNME
		* Ele vai criar um ambiente virtual e instalar as bibliotecas necessárias pra rodar o programa(sympy e matplotlib).

2. Executando o programa:
	No diretório do arquivo no terminal, rodar o comando:
		python3 metodos.py

3. Especificação da entrada:
	- A entrada é um arquivo txt que deve estar no mesmo diretório do projeto, com o nome entrada.txt.
	
	- Formato das entradas:
		- Há 3 formatos de entrada:
			(1) - nome_do_metodo y0 t0 h quantidade_de_passos função
			(2) - nome_do_metodo lista_ys_iniciais t0 h quantidade_de_passos função ordem
			(3) - nome_do_metodo y0 t0 h quantidade_de_passos função ordem
		
	- Códigos das entradas	
		* euler (1)
		* euler_inverso (1)
		* euler_aprimorado (1)
		* runge_kutta (1)
		* adam_bashforth (2)
		* adam_bashforth_by_euler (3)
		* adam_bashforth_by_euler_inverso (3)
		* adam_bashforth_by_euler_aprimorado (3)
		* adam_bashforth_by_runge_kutta (3)
		* adam_multon (2)
		* adam_multon_by_euler (3)
		* adam_multon_by_euler_inverso (3)
		* adam_multon_by_euler_aprimorado (3)
		* adam_multon_by_runge_kutta (3)
		* formula_inversa (2)
		* formula_inversa_by_euler (3)
		* formula_inversa_by_euler_inverso (3)
		* formula_inversa_by_euler_aprimorado (3)
		* formula_inversa_by_runge_kutta (3)
	

	- A entrada de y' deve ser especificando todos os sinais de operação:
		Ex: 1-t+4*y
4. Saída:
	- Um arquivo txt é gerado com as saidas do arquivo de entrada (saida.txt)
	
	
