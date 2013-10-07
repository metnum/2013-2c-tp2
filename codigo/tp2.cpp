// Escalante, Osinski, Raskovsky
// Grupo 16
#define grupo 16

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <iostream>
using namespace std;


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

void dibujar_diagonal(double * m, int n) {
    for (int i = 0; i < n; ++i) {
        printf("%f\n", pos(m, i, i));
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

    if (n == 2) {
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

        // Fill j2.y equation
        // f2.y + f3 + f5.y = 0
        pos(m, eq(2, "y"), f(2)) = y;
        pos(m, eq(2, "y"), f(3)) = 1;
        pos(m, eq(2, "y"), f(5)) = y;

        // Fill j(j_max).x equation
        // f(l_max).x + f(l_max - 1) = 0
        pos(m, eq(j_max, "x"), f(l_max)) = x;
        pos(m, eq(j_max, "x"), f(l_max - 1)) = 1;

        // Fill j(j_max).y equation
        // v1 = f(l_max).y
        pos(m, eq(j_max, "y"), v1_index) = 1;
        pos(m, eq(j_max, "y"), f(l_max)) = -y;
        return;
    }

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

void triangular_matriz(double * m, int n) {
    int i, j, max_i;
    double max_fabs, pivot, mult;
    // Itero por las filas
    for (i = 0; i < n; ++i) {
        // Busco el mayor pivot
        max_fabs = fabs(pos(m, i, i));
        max_i = i;
        for (j = i + 1; j < min(i + p + 1, n); ++j) {
            if (fabs(pos(m, j, i)) > max_fabs) {
                max_fabs = fabs(pos(m, j, i));
                max_i = j;
            }
        }
        // Si el mayor esta en otra fila permuto
        if (max_i != i) {
            permutar(m, i, max_i, n);
        }
        pivot = pos(m, i, i);
        for (j = i + 1; j < min(i + p + 1, n); ++j) {
            if (pos(m, j, i) != 0) {
                mult = - pos(m, j, i) / pivot;
                operacion(m, j, i, mult, n);
            }
        }
    }
}

void backwards_substitution(double * m, int n) {
    // Ir de atras para adelante, completando la diagonal
    double acumulado;
    int i, j;
    for (i = n -1; i >= 0; --i) {
        acumulado = posc(m, i);
        for (j = min(n - 1, i + p + q); j > i; --j) {
            acumulado -= pos(m, i, j) * pos(m, j, j);
        }
        pos(m, i, i) = acumulado / pos(m, i, i);
    }
}

double max_abs_diagonal(double * m, int n) {
    double max_abs = 0;
    for (int i = 0; i < n; ++i) {
        if (fabs(pos(m, i, i)) > max_abs) {
            max_abs = fabs(pos(m, i, i));
        }
    }
    return max_abs;
}

double procesar_seccion(double * m, double h, double span, double * C, int inicio, int fin) {

    int secciones = fin - inicio;

    // Inicializo con ceros, solo la parte de esta seccion
    for (int i = 0; i < banda * 4 * secciones; ++i) {
        m[i] = 0;
    }

    // Muevo el array C
    armar_matriz(m, secciones, h, span, &(C[inicio]));

    triangular_matriz(m, secciones * 4);
    backwards_substitution(m, secciones * 4);
    dibujar_diagonal(m, secciones * 4);
    return secciones;
}

double costo(int n, int secciones, double h, double span, double max_stress_seccion) {
    double x = span / n;
    double diagonal = sqrt( h * h + x * x );
    return ( diagonal * secciones + x * ( secciones * 2 - 2 ) + h * ( secciones - 1 ) ) * max_stress_seccion;
}

int main (int argc, char * argv[]) {

    // Todos los valores de entrada deben estar
    if (argc < 2 ) {
        cout << "Ejectuar ./tp2 <filename>" << endl;
        exit(-1);
    }

    // Leer de archivo y armar todo
    ifstream file (argv[1]);
    if (!file.is_open()) {
        cout << "No puedo abrir el archivo!!" << endl;
            return 1;
    }

    double span;
    double h;
    int n;
    double costo_pilar;
    double fmax;

    file >> span;
    file >> h;
    file >> n;

    // Agarrar memoria para la matriz banda y las cargas
    double * m = (double *) malloc(sizeof(double) * banda * 4 * n);  // 4n ecuaciones *  (banda + cargas Ci)
    double * C = (double *) malloc(sizeof(double) * (n - 1));  // vector de cargas iniciales

    if (m == NULL){
        cout << "ERROR: me quede sin memoria :( snif...\n" << endl;
        return 1;
    }

    // Lleno vector de cargas
    for(int i = 0; i < n - 1; i++) {
        file >> C[i];
    }

    file >> costo_pilar;
    file >> fmax;

    // Termine de leer del archivo

    // printf("span: %f, h: %f, n: %i, C: %f, fmax: %f\n", span, h, n, costo_pilar, fmax);


    // proceso la matriz
    int secciones = procesar_seccion(m, h, span, C, 0, 4);

    // obtengo el maximo stress
    double max_stress_seccion = max_abs_diagonal(m, secciones * 4);

    // calculo el costo
    double costo_seccion = costo(n, secciones, h, span, max_stress_seccion);

    cout << endl << max_stress_seccion << " " << costo_seccion << endl;

    // Empiezo a contar el tiempo
    init_time();
    double tiempo_total = get_time();

    // Cleanup
    file.close();
    free(m);
    free(C);
    return 0;
}

