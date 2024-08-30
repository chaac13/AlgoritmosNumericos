#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void gaussElimination(int dim, float** A, float* B);
void imprimirSolucoes(int dim, float* X, const char* nomeMetodo);

int main(int argc, char *argv[]) {
    int N_SISTEMAS, DIMENSAO;
    float PRECISAO;

    FILE *arquivo = fopen(argv[1], "r");
    if (arquivo == NULL) {
        printf("Erro ao abrir o arquivo.\n");
        return 1;
    }
    fscanf(arquivo, "%d %d %f\n", &N_SISTEMAS, &DIMENSAO, &PRECISAO);

    printf("N_SISTEMAS: %d\n", N_SISTEMAS);
    printf("DIMENSAO: %d\n", DIMENSAO);
    printf("PRECISAO: %f\n", PRECISAO);

    // Alocação dinâmica para a matriz A e o vetor B
    float** A = (float**)malloc(DIMENSAO * sizeof(float*));
    for (int i = 0; i < DIMENSAO; i++) {
        A[i] = (float*)malloc(DIMENSAO * sizeof(float));
    }
    float* B = (float*)malloc(DIMENSAO * sizeof(float));

    // Preenchendo a matriz A
    for (int i = 0; i < DIMENSAO; i++) {
        for (int j = 0; j < DIMENSAO; j++) {
            fscanf(arquivo, "%f", &A[i][j]);
        }
    }

    // Printando a matriz A
    printf("Matriz A:\n");
    for (int i = 0; i < DIMENSAO; i++) {
        for (int j = 0; j < DIMENSAO; j++) {
            printf("%.2f ", A[i][j]);
        }
        printf("\n");
    }

    for (int sistema = 0; sistema < N_SISTEMAS; sistema++) {
        // Lendo o vetor B
        for (int i = 0; i < DIMENSAO; i++) {
            fscanf(arquivo, "%f", &B[i]);
        }

        // Printando o vetor B
        printf("Vetor B:\n");
        for (int i = 0; i < DIMENSAO; i++) {
            printf("%.2f ", B[i]);
        }
        printf("\n");

        // Criando uma cópia da matriz A
        float** A_copia = (float**)malloc(DIMENSAO * sizeof(float*));
        for (int i = 0; i < DIMENSAO; i++) {
            A_copia[i] = (float*)malloc(DIMENSAO * sizeof(float));
            memcpy(A_copia[i], A[i], DIMENSAO * sizeof(float));  // Copia a linha inteira de A para A_copia
        }

        // Eliminação de Gauss
        clock_t start = clock();
        gaussElimination(DIMENSAO, A_copia, B);
        clock_t end = clock();
        double gaussTime = ((double)(end - start)) / CLOCKS_PER_SEC;
        printf("Tempo de execucao da Eliminacao de Gauss: %f segundos\n", gaussTime);

        // Liberando a memória da cópia da matriz A
        for (int i = 0; i < DIMENSAO; i++) {
            free(A_copia[i]);
        }
        free(A_copia);
    }

    // Liberando memória alocada
    for (int i = 0; i < DIMENSAO; i++) {
        free(A[i]);
    }
    free(A);
    free(B);

    fclose(arquivo);
    return 0;
}

void gaussElimination(int dim, float** A, float* B) {
    for (int i = 0; i < dim; i++) {
        for (int k = i + 1; k < dim; k++) {
            if (A[k][i] != 0) {
                float M = A[k][i] / A[i][i];
                for (int j = i; j < dim; j++) {
                    A[k][j] -= M * A[i][j];
                }
                B[k] -= M * B[i];
            }
        }
    }

    // Resolvendo o sistema após a eliminação (substituição retroativa)
    float* X = (float*)malloc(dim * sizeof(float));
    for (int i = dim - 1; i >= 0; i--) {
        X[i] = B[i];
        for (int j = i + 1; j < dim; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] /= A[i][i];
    }

    imprimirSolucoes(dim, X, "ELIMINACAO DE GAUSS");

    free(X);
}

void imprimirSolucoes(int dim, float* X, const char* nomeMetodo) {
    printf("\nSolucoes do sistema pelo metodo %s:\n", nomeMetodo);
    for (int i = 0; i < dim; i++) {
        printf("x%d = %.2f\n", i + 1, X[i]);
    }
}