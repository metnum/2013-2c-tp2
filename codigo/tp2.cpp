#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <float.h>
#include <math.h>
#include <sys/time.h>

#define p 4
#define q 4
#define banda (2 * p + q + 1 + 1)  // incluye la columna de cargas
#define pos_aux(A,i,j) A[(i) * banda + (j)]
#define pos(A, i, j) (pos_aux(A, i, j - i + p))
#define posc(A, i) (pos_aux(A, i, banda - 1))
#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) > (Y) ? (X) : (Y))

// auxiliares para armar la matriz
#define f(force) (force + 1)
#define eq(joint, axis) ((axis == "x") ? (2 * joint) : (2 * joint + 1))
#define h0_index  0
#define v0_index  1
#define v1_index  (4 * n - 1)

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
        printf("fila %2i.   ", i);
        for (int j = 0; j < n; ++j) {
            if (j - i + p >= 0 && j - i <= q + p) {
                printf("%05.2f  ", pos(m, i, j));
            } else {
                printf("00000  ");
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

void armar_matriz(double * m, int n, double h, double span, double * C) {

    double link_x = span / n;
    double link_y = h;
    double link_diag = pow(pow(link_x, 2) + pow(link_y, 2), 0.5);
    double x = link_y / link_diag;
    double y = link_x / link_diag;

    int k;  // iterator
    ////// Quantities
    int j_count = 2 * n;       // Number of joints
    int j_max = j_count - 1;   // Last joint index
    int l_count = 4 * n - 3;   // Number of links
    int l_max = l_count;    // Last link index

    // Fill j0.x equation
    // h0 = f1 + f2.x
    pos(m, eq(0, "x"), h0_index) = 1;
    pos(m, eq(0, "x"), f(1)) = -1;
    pos(m, eq(0, "x"), f(2)) = -x;

    // Fill j0.y equation
    // v0 = f2.y
    pos(m, eq(0, "y"), v0_index) = 1;
    pos(m, eq(0, "y"), f(2)) = -y;

    // Fill j1.x equation
    // f1 = f4
    pos(m, eq(1, "x"), f(1)) = 1;
    pos(m, eq(1, "x"), f(4)) = -1;

    // Fill j1.y equation
    // -c1 = f3
    pos(m, eq(1, "y"), f(3)) = 1;
    posc(m, eq(1, "y")) = C[0];

    // Fill j2.x equation
    // f2.x = f5.x + f6
    pos(m, eq(2, "x"), f(2)) = x;
    pos(m, eq(2, "x"), f(5)) = -x;
    pos(m, eq(2, "x"), f(6)) = -1;

    // Fill j2.y equation
    // f2.y + f3 + f5.y = 0;
    pos(m, eq(2, "y"), f(2)) = y;
    pos(m, eq(2, "y"), f(3)) = 1;
    pos(m, eq(2, "y"), f(5)) = y;

    // Fill j(n - 1).x: center lower.x equation
    // f(2 * n - 4) + f(2 * n - 3).x = f(2 * n) + f(2 * n + 1)
    pos(m, eq((n - 1), "x"), f(2 * n - 4)) = 1;
    pos(m, eq((n - 1), "x"), f(2 * n - 3)) = x;
    pos(m, eq((n - 1), "x"), f(2 * n)) = -1;
    pos(m, eq((n - 1), "x"), f(2 * n + 1)) = -x;

    // Fill j(n - 1).y: center lower.y equation
    // -c(n/2) = f(2 * n - 3).y + f(2 * n - 1) +  f(2 * n + 1).y
    posc(m, eq((n - 1), "y")) = C[n / 2 - 1];
    pos(m, eq((n - 1), "y"), f(2 * n - 3)) = y;
    pos(m, eq((n - 1), "y"), f(2 * n - 1)) = 1;
    pos(m, eq((n - 1), "y"), f(2 * n + 1)) = y;

    // Fill jn.x: center upper.x equation
    // f(2 * n - 2) = f(2 * n + 2)
    pos(m, eq(n, "x"), f(2 * n - 2)) = 1;
    pos(m, eq(n, "x"), f(2 * n + 2)) = -1;

    // Fill jn.y: center upper.y equation
    // f(2 * n - 1) = 0
    pos(m, eq(n, "y"), f(2 * n - 1)) = 1;

    // Fill j(j_max-2).x equation
    // f(l_max -5) = f(l_max - 1)
    pos(m, eq((j_max - 2), "x"), f(l_max - 5)) = 1;
    pos(m, eq((j_max - 2), "x"), f(l_max - 1)) = -1;

    // Fill j(j_max-2).y equation
    // -c(n-2) = f(l_max - 2)
    pos(m, eq((j_max - 2), "y"), f(l_max - 2)) = 1;
    posc(m, eq((j_max - 2), "y")) = C[n - 2];

    // Fill j(j_max-1).x equation
    // f(l_max -3)+ f(l_max - 4).x = f(l_max).x
    pos(m, eq((j_max - 1), "x"), f(l_max - 3)) = 1;
    pos(m, eq((j_max - 1), "x"), f(l_max - 4)) = x;
    pos(m, eq((j_max - 1), "x"), f(l_max)) = -x;

    // Fill j(j_max-1).y equation
    // f(l_max - 4).y + f(l_max - 2) + f(l_max).y = 0
    pos(m, eq((j_max - 1), "y"), f(l_max - 4)) = y;
    pos(m, eq((j_max - 1), "y"), f(l_max - 2)) = 1;
    pos(m, eq((j_max - 1), "y"), f(l_max)) = y;

    // Fill j(j_max).x equation
    // f(l_max).x + f(l_max - 1) = 0
    pos(m, eq(j_max, "x"), f(l_max)) = x;
    pos(m, eq(j_max, "x"), f(l_max - 1)) = 1;

    // Fill j(j_max).y equation
    // v1 = f(l_max).y
    pos(m, eq(j_max, "y"), v1_index) = 1;
    pos(m, eq(j_max, "y"), f(l_max)) = -y;

    // Fill j(2k) - 1, 1 < k < n / 2
    for (k = 3; k < n - 1; k += 2) {
        // f(2k -2) + f(2k -1).x = f(2k + 2)
        pos(m, eq(k, "x"), f(2 * k - 2)) = 1;
        pos(m, eq(k, "x"), f(2 * k - 1)) = x;
        pos(m, eq(k, "x"), f(2 * k + 2)) = -1;

        // -c(k/2) = f(2k - 1).y + f(2k + 1)
        posc(m, eq(k, "y")) = C[k / 2];
        pos(m, eq(k, "y"), f(2 * k - 1)) = y;
        pos(m, eq(k, "y"), f(2 * k + 1)) = 1;
    }

    // Fill j(2k)    , 1 < k < n / 2
    for (k = 4; k < n - 1; k += 2) {
        // f(2k - 2) = f(2k + 2) + f(2k + 1).x
        pos(m, eq(k, "x"), f(2 * k - 2)) = -1;
        pos(m, eq(k, "x"), f(2 * k + 2)) = 1;
        pos(m, eq(k, "x"), f(2 * k + 1)) = x;

        // f(2k - 1) + f(2k + 1).y = 0
        pos(m, eq(k, "y"), f(2 * k - 1)) = 1;
        pos(m, eq(k, "y"), f(2 * k + 1)) = y;
    }

    // Fill j(2k) - 1, n / 2 + 1 < k < n - 1
    for (k = n + 1; k < 2 * n - 3; k += 2) {
        // f(2k - 2) = f(2k + 2) + f(2k + 3).x
        pos(m, eq(k, "x"), f(2 * k - 2)) = -1;
        pos(m, eq(k, "x"), f(2 * k + 2)) = 1;
        pos(m, eq(k, "x"), f(2 * k + 3)) = x;

        // -c[k/2) = f(k * 2 + 1) + f(k * 2 + 3).y
        posc(m, eq(k, "y")) = C[k / 2];
        pos(m, eq(k, "y"), f(2 * k + 1)) = 1;
        pos(m, eq(k, "y"), f(2 * k + 3)) = y;
    }

    // Fill j(2k)    , n / 2 + 1 < k < n - 1
    for (k = n + 2; k < 2 * n - 3; k += 2) {
        // f(2k -3).x + f(2k - 2) = f(2k + 2)
        pos(m, eq(k, "x"), f(2 * k - 3)) = x;
        pos(m, eq(k, "x"), f(2 * k - 2)) = 1;
        pos(m, eq(k, "x"), f(2 * k + 2)) = -1;

        // f(2k -3).y + f(2k -1) = 0
        pos(m, eq(k, "y"), f(2 * k - 3)) = y;
        pos(m, eq(k, "y"), f(2 * k - 1)) = 1;
    }
    
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
    double * C = (double *) malloc(sizeof(double) * (n - 1));  // vector de cargas iniciales
    if (m == NULL){
        printf("ERROR: me quede sin memoria :( snif...\n");
        return 1;
    }

    // Inicializo con ceros
    for (int i = 0; i < banda * 4 * n; ++i) {
        m[i] = 0;
    }

    // Lleno la matriz original con algo
    // for (int i = 0; i < n * 4; ++i) {
    //     for (int j = 0; j < n * 4; ++j) {
    //         if (j - i + p >= 0 && j - i <= q) {
    //             if (i <= j)
    //             pos(m, i, j) = i + j + 1;
    //         }
    //     }
    //     // y los C
    //     posc(m, i) = i;
    // }

    // Seteo la diagonal
    // for (int i = 0; i < 4 * n; ++i) {
    //     for (int j = 0; j < 4 * n; ++j) {
    //         if (i == j) {
    //             pos(m, i, j) = 8;
    //         }
    //     }
    // }

    // Obtengo el vector de cargas
    for (int i = 0; i < n - 1; ++i) {
        C[i] = i + 1;
        printf("%f ", C[i]);
    }
    printf("\n");

    armar_matriz(m, n, h, span, C);

    dibujar_matriz(m, 4 * n);
    printf("\n");
    // permutar(m, 5, 7, 4 * n);
    // operacion(m, 11, 9, 1, 4 * n);
    // dibujar_matriz(m, 4 * n);


    // Empiezo a contar el tiempo
    init_time();
    double tiempo_total = get_time();

    free(m);
    return 0;
}

