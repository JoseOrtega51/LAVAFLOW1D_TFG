# LAVAFLOW_1D
 Código para TFG de física: Predicción del avance de coladas de lava mediante simulación numérica
 
- LAVAFLOW_1D.c: Archivo principal de simulación.
- Carpeta lib: Contiene un fichero de funciones y otro de configuración para complementar LAVAFLOW_1D.c
- FIT_BINGHAM_PYTHON: Carpeta que incluye los archivos para calcular y ajustar el modelo reológico para un determinado perefil de temperatura.
- gif.py: Realiza gifs de la simulación.

 Puede usar los siguientes modelos reológicos:
  * Fluido Newtoniano viscosidad constante. (1)
  * Fluido Bingham con viscosidad y esfuerzo umbral constantes. (2)
  * Fluido Bingham con viscosidad y esfuerzo umbral variables asumiendo una temperatura constante en la columna. (3)
  * Fluido Bingham con viscosidad y esfuerzo umbral variables asumiendo un perfil de temperatura ajustable. (4)
