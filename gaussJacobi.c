
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void gaussJacobi(int dim, float** A, float* B, float* X, float precisao);
void lerDados(FILE *arquivo, int *N_SISTEMAS, int *DIMENSAO, float *PRECISAO, float*** A, float** B);
void imprimirSolucoes(int dim, float* X, const char* nomeMetodo);

int main(int argc, char *argv[]) {
    int N_SISTEMAS, DIMENSAO;
    float PRECISAO;
    float **A = NULL;
    float *B = NULL;

    FILE *arquivo = fopen(argv[1], "r");
    if (arquivo == NULL) {
        printf("Erro ao abrir o arquivo.\n");
        return 1;
    }

    lerDados(arquivo, &N_SISTEMAS, &DIMENSAO, &PRECISAO, &A, &B);

    for (int sistema = 0; sistema < N_SISTEMAS; sistema++) {
        // Vetor X para armazenar as soluções
        float *X = (float*)calloc(DIMENSAO, sizeof(float));
        if (X == NULL) {
            printf("Erro ao alocar memória para o vetor X.\n");
            fclose(arquivo);
            return 1;
        }

        // Lendo o vetor B para o sistema atual
        for (int i = 0; i < DIMENSAO; i++) {
            fscanf(arquivo, "%f", &B[i]);
        }

        clock_t start = clock();
        gaussJacobi(DIMENSAO, A, B, X, PRECISAO);
        clock_t end = clock();

        imprimirSolucoes(DIMENSAO, X, "GAUSS-JACOBI");
        printf("Tempo de execucao do Metodo de Gauss-Jacobi: %f segundos\n", ((double)(end - start)) / CLOCKS_PER_SEC);

        free(X);
    }

    for (int i = 0; i < DIMENSAO; i++) {
        free(A[i]);
    }
    free(A);
    free(B);
    fclose(arquivo);

    return 0;
}

void gaussJacobi(int dim, float** A, float* B, float* X, float precisao) {
    float *x_1 = (float*)calloc(dim, sizeof(float));  // x(k)
    float *x_2 = (float*)calloc(dim, sizeof(float));  // x(k+1)

    int PRIMEIRA_ITERACAO = 1;
    float maior_dr = 9999;

    while (maior_dr > precisao) {
        if (PRIMEIRA_ITERACAO) {
            for (int i = 0; i < dim; i++) {
                x_1[i] = B[i] / A[i][i];
            }
            PRIMEIRA_ITERACAO = 0;
        } else {
            for (int i = 0; i < dim; i++) {
                x_1[i] = x_2[i];
            }
        }

        for (int i = 0; i < dim; i++) {
            x_2[i] = B[i];
            for (int j = 0; j < dim; j++) {
                if (i != j) {
                    x_2[i] -= A[i][j] * x_1[j];
                }
            }
            x_2[i] /= A[i][i];
        }

        float maior_d = 0;
        for (int i = 0; i < dim; i++) {
            if (fabs(x_2[i] - x_1[i]) > maior_d) {
                maior_d = fabs(x_2[i] - x_1[i]);
            }
        }

        float maior_x_2 = 0;
        for (int i = 0; i < dim; i++) {
            if (x_2[i] > maior_x_2) {
                maior_x_2 = x_2[i];
            }
        }

        maior_dr = maior_d / maior_x_2;
    }

    for (int i = 0; i < dim; i++) {
        X[i] = x_2[i];
    }
    free(x_1);
    free(x_2);
}

void lerDados(FILE *arquivo, int *N_SISTEMAS, int *DIMENSAO, float *PRECISAO, float*** A, float** B) {
    fscanf(arquivo, "%d %d %f", N_SISTEMAS, DIMENSAO, PRECISAO);

    *A = (float**)malloc((*DIMENSAO) * sizeof(float*));
    for (int i = 0; i < *DIMENSAO; i++) {
        (*A)[i] = (float*)malloc((*DIMENSAO) * sizeof(float));
    }

    *B = (float*)malloc((*DIMENSAO) * sizeof(float));

    for (int i = 0; i < *DIMENSAO; i++) {
        for (int j = 0; j < *DIMENSAO; j++) {
            fscanf(arquivo, "%f", &(*A)[i][j]);
        }
    }
}
void imprimirSolucoes(int dim, float* X, const char* nomeMetodo) {
    printf("Solucoes do sistema pelo metodo %s:\n", nomeMetodo);
    for (int i = 0; i < dim; i++) {
        printf("x%d = %.6f\n", i + 1, X[i]);
    }
    printf("\n");
}