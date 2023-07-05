Contenido de la carpeta:
- compute_u_mean.py: Calcula las velocidades promedio a partir de un perfil de temperatura y los parámetros para las leyes reológicas.
- fit_u_mean.py: Ajusta los las funciones a partir de los datos de compute_u_mean.py y saca un fichero de texto DATA_FIT_BINGHAM.txt con los resultados.
- Diffusion.c: Calcular numéricamente los perfiles de temperatura con el modelo de difusión de temperatura.
- Ts.py: Ajusta la temperatura en la superficie en función de la temperatura media para el modelo de difusión.
- Carpeta lib: Tiene ficheros con funciones.
  
El proceso para ajustar los parámetros de una reología con perfil de temperatura es el siguiente:
1. (Si se va a emplear un modelo difusivo). Diffusion.c 
1. compute_u_mean.py
1. fit_u_mean.py
1. Copiar DATA_FIT_BINGHAM.txt en la carpeta que contiene a LAVAFLOW_1D.c
