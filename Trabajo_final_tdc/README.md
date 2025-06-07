Por Juan Diego Gómez Guzmán
# Simulador de Dinámica Molecular con Lennard-Jones (Paralelo)

Este repositorio contiene el código fuente para un simulador de dinámica molecular en 3D que modela el comportamiento de partículas (como átomos de Argón) interactuando bajo el potencial de Lennard-Jones. El proyecto incluye tanto una versión secuencial de referencia como una versión paralela optimizada con OpenMP para aprovechar las arquitecturas de múltiples núcleos.

## Características

* **Modelo Físico:** Implementa el potencial de **Lennard-Jones** para calcular las fuerzas de repulsión y atracción entre partículas.
* **Dinámica:** Utiliza el integrador numérico **Velocity Verlet** para una evolución temporal estable y con buena conservación de la energía.
* **Geometría:** Emplea **Condiciones de Frontera Periódicas (PBC)** con la convención de la imagen mínima para simular un sistema macroscópico.
* **Flexibilidad:** La simulación es completamente configurable a través de archivos de texto. Se pueden modificar el número de partículas, tamaño de la caja, duración, propiedades de la materia, etc., sin necesidad de recompilar.
* **Paralelismo:** La versión paralela utiliza **OpenMP** para acelerar significativamente el cálculo de fuerzas, que es el cuello de botella computacional del sistema.

---

## Prerrequisitos

Para compilar y ejecutar este proyecto, necesitas un compilador de C/C++ que soporte OpenMP. La configuración recomendada es:
* **Sistema Operativo:** Windows, Linux o macOS.
* **Compilador:** `gcc` (parte de GCC o MinGW-w64 en Windows).
* **Herramienta de construcción:** `make`.

En Windows, se recomienda encarecidamente instalar **MSYS2** y usar la terminal **UCRT64** para asegurar que `gcc` y `make` estén disponibles y configurados correctamente.

---

## Compilación

El proyecto utiliza un `Makefile` que automatiza todo el proceso de compilación.

1.  Abre una terminal en la **carpeta raíz** del proyecto (donde se encuentra el `Makefile`).
2.  Ejecuta uno de los siguientes comandos:

    * **Para compilar ambas versiones (secuencial y paralela):**
        ```bash
        make
        ```
        o
        ```bash
        make all
        ```

    * **Para compilar únicamente la versión paralela:**
        ```bash
        make simulador_paralelo
        ```

    * **Para limpiar todos los archivos compilados:**
        ```bash
        make clean
        ```

Los archivos ejecutables se crearán dentro de la carpeta `bin/`.

---

## Uso

Una vez compilado, el programa se ejecuta desde la terminal, especificando el archivo de configuración del escenario que se desea simular.

#### Ejecución de un Escenario

* **Versión Secuencial:**
    ```bash
    ./bin/simulador_secuencial config/benchmark.txt
    ```

* **Versión Paralela:**
    ```bash
    ./bin/simulador_paralelo config/benchmark.txt
    ```

#### Pruebas de Rendimiento y Escalabilidad

Para medir el *speed-up* y la eficiencia, puedes controlar el número de hilos que utiliza la versión paralela mediante la variable de entorno `OMP_NUM_THREADS`.

**Ejemplo en terminales bash (Linux, macOS, MSYS2):**
```bash
# Ejecutar con 1 hilo (base para T1)
export OMP_NUM_THREADS=1
./bin/simulador_paralelo config/benchmark_N512.txt

# Ejecutar con 4 hilos
export OMP_NUM_THREADS=4
./bin/simulador_paralelo config/benchmark_N512.txt

# Ejecutar con 8 hilos
export OMP_NUM_THREADS=8
./bin/simulador_paralelo config/benchmark_N512.txt

Ejemplo en PowerShell (Windows):

$env:OMP_NUM_THREADS=4
.\bin\simulador_paralelo.exe config\benchmark_N512.txt

Estructura del Repositorio
/
├── bin/          # Contiene los programas ejecutables compilados (ignorado por Git).
├── config/       # Almacena los archivos de configuración para diferentes escenarios.
├── doc/          # Documentación del proyecto (informes, etc.).
├── results/      # Salidas de las simulaciones (logs, datos) (ignorado por Git).
├── src/          # Contiene todo el código fuente (.c).
├── .gitignore    # Especifica los archivos y carpetas que Git debe ignorar.
├── Makefile      # La "receta" para compilar el proyecto.
└── README.md     # Este archivo.
