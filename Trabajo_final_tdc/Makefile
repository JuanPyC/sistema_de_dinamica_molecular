# Makefile para el Proyecto de Dinámica Molecular

# --- Variables de Compilación ---
# Compilador que usaremos
CC = gcc

# Banderas para la compilación. -O2 es optimización, -Wall muestra todas las advertencias.
# -fopenmp es para activar OpenMP.
CFLAGS = -O2 -Wall -fopenmp

# Banderas para el enlazador (linker). Le decimos que use las librerías de OpenMP y la de matemáticas (-lm).
# Añadimos la bandera robusta para asegurar que sea una aplicación de consola en Windows.
LDFLAGS = -fopenmp -lm -Wl,-subsystem,console

# --- Variables de Directorios y Archivos ---
# Directorio del código fuente
SRCDIR = src
# Directorio de los ejecutables compilados
BINDIR = bin

# Nombres de los archivos fuente
SRC_SEQ = $(SRCDIR)/simulador_refactorizado.c
SRC_PAR = $(SRCDIR)/simulador_paralelo.c

# Nombres de los ejecutables de destino
TARGET_SEQ = $(BINDIR)/simulador_secuencial
TARGET_PAR = $(BINDIR)/simulador_paralelo


# --- Reglas de Construcción (Las "Recetas") ---

# La regla 'all' es la que se ejecuta por defecto si solo escribes "make".
# Le decimos que construya tanto la versión secuencial como la paralela.
all: $(TARGET_SEQ) $(TARGET_PAR)

# Regla para construir el ejecutable PARALELO
$(TARGET_PAR): $(SRC_PAR)
	@echo "Compilando la versión paralela..."
	@mkdir -p $(BINDIR) # Crea la carpeta bin/ si no existe
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)
	@echo "¡Versión paralela compilada en $(TARGET_PAR)!"

# Regla para construir el ejecutable SECUENCIAL
# Nota: La única diferencia real es el archivo fuente que usa.
$(TARGET_SEQ): $(SRC_SEQ)
	@echo "Compilando la versión secuencial..."
	@mkdir -p $(BINDIR) # Crea la carpeta bin/ si no existe
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)
	@echo "¡Versión secuencial compilada en $(TARGET_SEQ)!"

# Regla 'clean' para limpiar los archivos compilados
clean:
	@echo "Limpiando el directorio de compilación..."
	@rm -f $(BINDIR)/*
	@echo "Limpieza completada."

# Le dice a 'make' que 'all' y 'clean' no son archivos, sino nombres de comandos.
.PHONY: all clean