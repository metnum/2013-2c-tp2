#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <float.h>
#include <math.h>
#include <sys/time.h>

#define p 2
#define q 2
#define banda (2 * p + q + 1 + 1)  // incluye la columna de cargas
#define pos_aux(A,i,j) A[(i) * banda + (j)]
#define pos(A, i, j) (pos_aux(A, i, j - i + p))
#define posc(A, i) (pos_aux(A, i, banda - 1))
#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) > (Y) ? (X) : (Y))

// Funciones para calcular el tiempo
struct timeval start, end;
void init_time()
{
    gettimeofday(&start,NULL);
}

double get_time()
{
    gettimeofday(&end,NULL);
    return (1000000*(end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec))/1000000.0;
}

void dibujar_matriz(double * m, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j - i + p >= 0 && j - i <= q + p) {
                printf("%03.0f ", pos(m, i, j));
            } else {
                printf("000 ");
            }
        }
        printf("%06.2f\n", posc(m, i));
    }
}

void permutar(double * m, int i, int j, int n) {
    // i must be lower than j
    double swap;
    for (int k = i; k < min(n, i + p + q + 1); ++k) {
        if (k < i + q + 1) {
            swap = pos(m, i, k);
        } else {
            swap = 0;
        }
        pos(m, i, k) = pos(m, j, k);
        pos(m, j, k) = swap;
    } 
    swap = posc(m, i);
    posc(m, i) = posc(m, j);
    posc(m, j) = swap;
}

void operacion(double * m, int i, int j, double mult, int n) {
    // j must be lower than i
    //m[i] = m[i] + m[j] * mult
    for (int k = max(j - p, 0); k < min(n, j + p + q + 1); ++k) {
        printf("poniendo en la fila i %i, pos k %i, fila j %i * mult %f\n", i, k, j, mult);
        pos(m, i, k) = pos(m, i, k) + mult * pos(m, j, k);
    } 
    posc(m, i) = posc(m, i) + mult * posc(m, j);
}

int main (int argc, char * argv[]) {

    // Todos los valores de entrada deben estar
    if (argc < 4 ) {
        printf("Ejectuar ./tp2 <span> <h> <n>\n");
        exit(-1);
    }

    // Leo la entrada
    char * end;
    double span = strtod(argv[1], &end);
    if (argv[1] == end) {
        printf("Error: El <span> debe ser un double\n");
        exit(-1);
    }
    double h = strtod(argv[2], &end);
    if (argv[2] == end) {
        printf("Error: El <h> debe ser un double\n");
        exit(-1);
    }
    double n = strtod(argv[3], &end);
    if (argv[3] == end) {
        printf("Error: El <n> debe ser un double\n");
        exit(-1);
    }

    printf("span: %f, h: %f, n: %f\n", span, h, n);

    double * m = (double *) malloc(sizeof(double) * banda * 4 * n);  // 4n ecuaciones *  (banda + cargas Ci)
    if (m == NULL){
        printf("ERROR: me quede sin memoria :( snif...\n");
        return 1;
    }

    // Inicializo con ceros
    for (int i = 0; i < banda * 4 * n; ++i) {
        m[i] = 0;
    }

    // Lleno la matriz original con algo
    for (int i = 0; i < n * 4; ++i) {
        for (int j = 0; j < n * 4; ++j) {
            if (j - i + p >= 0 && j - i <= q) {
                pos(m, i, j) = i + 1;
            }
        }
        // y los C
        posc(m, i) = i;
    }

    // Seteo la diagonal
    // for (int i = 0; i < 4 * n; ++i) {
    //     for (int j = 0; j < 4 * n; ++j) {
    //         if (i == j) {
    //             pos(m, i, j) = 8;
    //         }
    //     }
    // }

    dibujar_matriz(m, 4 * n);
    printf("\n");
    //permutar(m, 5, 7, 4 * n);
    operacion(m, 7, 5, -1.5, 4 * n);
    dibujar_matriz(m, 4 * n);


    // Empiezo a contar el tiempo
    init_time();
    double tiempo_total = get_time();

    free(m);
    return 0;
}

