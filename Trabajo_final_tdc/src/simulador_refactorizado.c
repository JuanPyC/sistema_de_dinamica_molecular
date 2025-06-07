#include <stdio.h>    // Para funciones de entrada/salida estándar (printf, fopen, etc.)
#include <stdlib.h>   // Para manejo de memoria dinámica (malloc, free, exit)
#include <math.h>     // Para funciones matemáticas (sqrt, pow, cbrt)
#include <string.h>   // Para manipulación de cadenas de texto (strcmp)
#include <omp.h>      // Necesario para las funciones y directivas de OpenMP (paralelismo)

// --- Estructuras de Datos ---

// Define un tipo de dato 'Vector' para manejar coordenadas 3D de forma sencilla.
typedef struct { double x, y, z; } Vector;

// Define un tipo de dato 'Particle' que contiene toda la información de una partícula.
typedef struct { 
    int id;             // Identificador único para cada partícula.
    Vector pos;         // Vector de posición (x, y, z).
    Vector vel;         // Vector de velocidad (vx, vy, vz).
    Vector force;       // Vector de la fuerza neta que actúa sobre la partícula.
} Particle;

// Define una estructura 'Config' para agrupar todos los parámetros de la simulación.
// Esto facilita pasar la configuración completa a las funciones.
typedef struct {
    int N_PARTICLES, N_STEPS, PRINT_FREQ;
    double BOX_SIZE, TIMESTEP, MASS, EPSILON, SIGMA, CUTOFF_RADIUS;
} Config;

// --- Prototipos de Funciones ---
// Declaraciones anticipadas de todas las funciones que se usarán en el programa.
void load_config(Config *config, const char* filename);
void initialize_particles_lattice(Particle particles[], const Config* config);
double calculate_forces_and_potential_parallel(Particle particles[], const Config* config);
void update_positions(Particle particles[], const Config* config);
void update_velocities(Particle particles[], const Vector old_forces[], const Config* config);
void apply_periodic_boundaries(Particle particles[], const Config* config);
double calculate_kinetic_energy(Particle particles[], const Config* config);


int main(int argc, char *argv[]) {
    // Valida que el programa se haya ejecutado con exactamente un argumento (el archivo de config).
    if (argc != 2) {
        fprintf(stderr, "Uso: %s <archivo_config>\n", argv[0]);
        return 1; // Termina con un código de error si el uso es incorrecto.
    }
    
    // Almacena la configuración de la simulación.
    Config config;
    // Carga los parámetros desde el archivo especificado en la línea de comandos.
    load_config(&config, argv[1]);

    // Reserva memoria en el heap para el array de partículas, basado en N_PARTICLES.
    Particle *particles = (Particle*) malloc(config.N_PARTICLES * sizeof(Particle));
    // Reserva memoria para guardar las fuerzas del paso anterior, necesario para Verlet.
    Vector *old_forces = (Vector*) malloc(config.N_PARTICLES * sizeof(Vector));
    
    // Verifica si la reserva de memoria falló.
    if (!particles || !old_forces) {
        fprintf(stderr, "Error: Malloc falló para N=%d.\n", config.N_PARTICLES); 
        return 1;
    }

    // Coloca las partículas en una red cúbica inicial.
    initialize_particles_lattice(particles, &config);

    // --- Inicio de la Simulación ---
    printf("Iniciando simulación PARALELA con archivo: %s\n", argv[1]);
    // Imprime el número máximo de hilos que OpenMP puede usar (usualmente, el número de núcleos de la CPU).
    printf("Usando un máximo de %d hilos (threads).\n", omp_get_max_threads());
    printf("N=%d, Pasos=%d, dt=%.4f\n", config.N_PARTICLES, config.N_STEPS, config.TIMESTEP);
    printf("------------------------------------------------------------------------\n");
    printf("Paso\tT_Total(s)\tE_Cinetica\tE_Potencial\tE_Total\n");
    printf("------------------------------------------------------------------------\n");

    // Inicia el cronómetro.
    double start_time = omp_get_wtime();
    // Calcula las fuerzas y la energía potencial iniciales (en t=0).
    double potential_energy = calculate_forces_and_potential_parallel(particles, &config);

    // Bucle principal de la simulación que se ejecuta N_STEPS veces.
    for (int step = 1; step <= config.N_STEPS; step++) {
        
        // --- INICIO DEL ALGORITMO VELOCITY VERLET ---
        
        // 1. Actualiza las posiciones r(t+dt) usando las fuerzas actuales F(t).
        update_positions(particles, &config);
        
        // 2. Guarda las fuerzas actuales F(t) en 'old_forces' antes de que se sobreescriban.
        for(int i=0; i<config.N_PARTICLES; i++) old_forces[i] = particles[i].force;

        // 3. Aplica condiciones de frontera a las nuevas posiciones para mantenerlas dentro de la caja.
        apply_periodic_boundaries(particles, &config);

        // 4. Calcula las nuevas fuerzas F(t+dt) basadas en las nuevas posiciones r(t+dt).
        potential_energy = calculate_forces_and_potential_parallel(particles, &config);

        // 5. Actualiza las velocidades v(t+dt) usando el promedio de las fuerzas vieja y nueva.
        update_velocities(particles, old_forces, &config);

        // --- FIN DEL ALGORITMO VELOCITY VERLET ---

        // Imprime el estado del sistema cada 'PRINT_FREQ' pasos.
        if (step % config.PRINT_FREQ == 0) {
            double kinetic_energy = calculate_kinetic_energy(particles, &config);
            double total_energy = kinetic_energy + potential_energy;
            double elapsed_time = omp_get_wtime() - start_time;
            // Imprime los datos en columnas para un fácil análisis.
            printf("%d\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n", step, elapsed_time, kinetic_energy, potential_energy, total_energy);
        }
    }

    // Detiene el cronómetro e imprime el tiempo total de ejecución.
    double end_time = omp_get_wtime();
    printf("------------------------------------------------------------------------\n");
    printf("Simulación completada en %.4f segundos.\n", end_time - start_time);

    // Libera la memoria reservada para evitar fugas de memoria.
    free(particles);
    free(old_forces);
    return 0;
}

// Carga los parámetros desde un archivo de texto y los almacena en la estructura Config.
void load_config(Config *config, const char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error: No se pudo abrir el archivo de configuración %s\n", filename);
        exit(1);
    }
    char line[128], key[64];
    double value;
    // Lee el archivo línea por línea.
    while (fgets(line, sizeof(line), file)) {
        // Ignora líneas vacías o que comiencen con '#'.
        if (line[0] == '#' || line[0] == '\n') continue;
        // Parsea la línea para obtener una clave (string) y un valor (double).
        if (sscanf(line, "%s %lf", key, &value) == 2) {
            // Compara la clave y asigna el valor al campo correspondiente de la estructura.
            if (strcmp(key, "N_PARTICLES") == 0) config->N_PARTICLES = (int)value;
            else if (strcmp(key, "N_STEPS") == 0) config->N_STEPS = (int)value;
            else if (strcmp(key, "PRINT_FREQ") == 0) config->PRINT_FREQ = (int)value;
            else if (strcmp(key, "BOX_SIZE") == 0) config->BOX_SIZE = value;
            else if (strcmp(key, "TIMESTEP") == 0) config->TIMESTEP = value;
            else if (strcmp(key, "MASS") == 0) config->MASS = value;
            else if (strcmp(key, "EPSILON") == 0) config->EPSILON = value;
            else if (strcmp(key, "SIGMA") == 0) config->SIGMA = value;
        }
    }
    // Calcula un parámetro derivado (radio de corte) basado en sigma.
    config->CUTOFF_RADIUS = 2.5 * config->SIGMA;
    fclose(file); // Cierra el archivo.
}

// Inicializa las posiciones de las partículas en una red cúbica simple.
void initialize_particles_lattice(Particle particles[], const Config* config) {
    int n_per_dim = (int)cbrt(config->N_PARTICLES); // Calcula cuántas partículas hay por dimensión.
    double spacing = config->BOX_SIZE / n_per_dim; // Calcula el espaciado entre ellas.
    int p_idx = 0;
    for (int x = 0; x < n_per_dim; x++) {
        for (int y = 0; y < n_per_dim; y++) {
            for (int z = 0; z < n_per_dim; z++) {
                if (p_idx < config->N_PARTICLES) {
                    particles[p_idx].id = p_idx;
                    // Asigna la posición en el centro de cada celda de la red.
                    particles[p_idx].pos = (Vector){(x + 0.5) * spacing, (y + 0.5) * spacing, (z + 0.5) * spacing};
                    // Inicializa velocidades y fuerzas a cero.
                    particles[p_idx].vel = (Vector){0.0, 0.0, 0.0};
                    particles[p_idx].force = (Vector){0.0, 0.0, 0.0};
                    p_idx++;
                }
            }
        }
    }
}

// Calcula las fuerzas y la energía potencial del sistema de forma paralela.
// Este es el "cuello de botella" computacional del programa.
double calculate_forces_and_potential_parallel(Particle particles[], const Config* config) {
    // Pre-calcula el radio de corte al cuadrado para evitar usar sqrt en el bucle.
    double r_cut_sq = config->CUTOFF_RADIUS * config->CUTOFF_RADIUS;
    double potential_energy = 0.0; // Variable para acumular la energía potencial total.

    // Primero, reinicia las fuerzas de todas las partículas a cero. (Este bucle es rápido y no se paraleliza).
    for (int i = 0; i < config->N_PARTICLES; i++) {
        particles[i].force = (Vector){0.0, 0.0, 0.0};
    }

    // Directiva OpenMP: Inicia una región paralela y distribuye las iteraciones del bucle 'i'.
    // 'reduction(+:potential_energy)': Soluciona la condición de carrera para la energía potencial.
    // Cada hilo tendrá su propia copia de 'potential_energy' y al final se sumarán todas.
    #pragma omp parallel for reduction(+:potential_energy)
    for (int i = 0; i < config->N_PARTICLES; i++) {
        // El bucle interno itera sobre los pares (i, j) sin repetir (j = i + 1).
        for (int j = i + 1; j < config->N_PARTICLES; j++) {
            Vector dr; // Vector de distancia entre la partícula i y j.
            dr.x = particles[i].pos.x - particles[j].pos.x;
            dr.y = particles[i].pos.y - particles[j].pos.y;
            dr.z = particles[i].pos.z - particles[j].pos.z;

            // --- Convención de Imagen Mínima (para Condiciones de Frontera Periódicas) ---
            // Si la distancia en una dimensión es mayor que la mitad de la caja, se corrige
            // para usar la distancia a la "imagen" de la partícula en la caja vecina.
            if (dr.x > config->BOX_SIZE/2.0) dr.x -= config->BOX_SIZE;
            else if (dr.x < -config->BOX_SIZE/2.0) dr.x += config->BOX_SIZE;
            if (dr.y > config->BOX_SIZE/2.0) dr.y -= config->BOX_SIZE;
            else if (dr.y < -config->BOX_SIZE/2.0) dr.y += config->BOX_SIZE;
            if (dr.z > config->BOX_SIZE/2.0) dr.z -= config->BOX_SIZE;
            else if (dr.z < -config->BOX_SIZE/2.0) dr.z += config->BOX_SIZE;
            
            // Calcula la distancia al cuadrado.
            double r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

            // Solo calcula la interacción si el par está dentro del radio de corte.
            if (r2 < r_cut_sq) {
                // Cálculo eficiente de los términos de Lennard-Jones.
                double r6_inv = pow(config->SIGMA*config->SIGMA/r2, 3);
                double r12_inv = r6_inv * r6_inv;
                
                // Acumula la energía potencial de este par.
                potential_energy += 4.0 * config->EPSILON * (r12_inv - r6_inv);
                
                // Calcula la magnitud de la fuerza dividida por la distancia (F/r).
                double force_mag = (24.0 * config->EPSILON / r2) * (2.0 * r12_inv - r6_inv);
                
                // Calcula el vector de fuerza.
                Vector force_vec = {force_mag * dr.x, force_mag * dr.y, force_mag * dr.z};
                
                // --- Actualización Atómica de Fuerzas ---
                // Se usa '#pragma omp atomic' para evitar condiciones de carrera.
                // Asegura que la operación de suma se complete sin que otro hilo interfiera,
                // lo cual es crucial cuando dos hilos intentan actualizar la misma partícula 'j'.
                #pragma omp atomic
                particles[i].force.x += force_vec.x;
                #pragma omp atomic
                particles[i].force.y += force_vec.y;
                #pragma omp atomic
                particles[i].force.z += force_vec.z;
                
                // Por la Tercera Ley de Newton, la fuerza en j es opuesta a la de i.
                #pragma omp atomic
                particles[j].force.x -= force_vec.x;
                #pragma omp atomic
                particles[j].force.y -= force_vec.y;
                #pragma omp atomic
                particles[j].force.z -= force_vec.z;
            }
        }
    }
    return potential_energy; // Devuelve la energía potencial total del sistema.
}

// Actualiza las posiciones de las partículas usando la primera parte de la fórmula de Verlet.
void update_positions(Particle particles[], const Config* config) {
    for (int i = 0; i < config->N_PARTICLES; i++) {
        // Calcula la aceleración actual (a = F/m).
        double accel_x = particles[i].force.x / config->MASS;
        double accel_y = particles[i].force.y / config->MASS;
        double accel_z = particles[i].force.z / config->MASS;

        // r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
        particles[i].pos.x += particles[i].vel.x * config->TIMESTEP + 0.5 * accel_x * config->TIMESTEP * config->TIMESTEP;
        particles[i].pos.y += particles[i].vel.y * config->TIMESTEP + 0.5 * accel_y * config->TIMESTEP * config->TIMESTEP;
        particles[i].pos.z += particles[i].vel.z * config->TIMESTEP + 0.5 * accel_z * config->TIMESTEP * config->TIMESTEP;
    }
}

// Actualiza las velocidades usando la segunda parte de la fórmula de Verlet.
void update_velocities(Particle particles[], const Vector old_forces[], const Config* config) {
    for (int i = 0; i < config->N_PARTICLES; i++) {
        // Calcula la aceleración del paso anterior (a(t)).
        double old_accel_x = old_forces[i].x / config->MASS;
        double old_accel_y = old_forces[i].y / config->MASS;
        double old_accel_z = old_forces[i].z / config->MASS;
        
        // Calcula la nueva aceleración (a(t+dt)) a partir de las fuerzas ya actualizadas.
        double new_accel_x = particles[i].force.x / config->MASS;
        double new_accel_y = particles[i].force.y / config->MASS;
        double new_accel_z = particles[i].force.z / config->MASS;

        // v(t+dt) = v(t) + 0.5 * (a(t) + a(t+dt)) * dt
        particles[i].vel.x += 0.5 * (old_accel_x + new_accel_x) * config->TIMESTEP;
        particles[i].vel.y += 0.5 * (old_accel_y + new_accel_y) * config->TIMESTEP;
        particles[i].vel.z += 0.5 * (old_accel_z + new_accel_z) * config->TIMESTEP;
    }
}

// Aplica las condiciones de frontera periódicas.
void apply_periodic_boundaries(Particle particles[], const Config* config) {
    for (int i = 0; i < config->N_PARTICLES; i++) {
        // Si la partícula se sale por un lado de la caja, reaparece por el lado opuesto.
        if (particles[i].pos.x >= config->BOX_SIZE) particles[i].pos.x -= config->BOX_SIZE;
        else if (particles[i].pos.x < 0) particles[i].pos.x += config->BOX_SIZE;
        
        if (particles[i].pos.y >= config->BOX_SIZE) particles[i].pos.y -= config->BOX_SIZE;
        else if (particles[i].pos.y < 0) particles[i].pos.y += config->BOX_SIZE;

        if (particles[i].pos.z >= config->BOX_SIZE) particles[i].pos.z -= config->BOX_SIZE;
        else if (particles[i].pos.z < 0) particles[i].pos.z += config->BOX_SIZE;
    }
}

// Calcula la energía cinética total del sistema.
double calculate_kinetic_energy(Particle particles[], const Config* config) {
    double kinetic_energy = 0.0;
    // Suma la energía cinética (0.5 * m * v^2) de cada partícula.
    for (int i = 0; i < config->N_PARTICLES; i++) {
        double v2 = particles[i].vel.x * particles[i].vel.x +
                      particles[i].vel.y * particles[i].vel.y +
                      particles[i].vel.z * particles[i].vel.z;
        kinetic_energy += 0.5 * config->MASS * v2;
    }
    return kinetic_energy;
}