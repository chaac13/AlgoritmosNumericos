#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void luDecomposition(int dim, float** A, float* B);
void forwardSubstitution(int dim, float** L, float* B, float* y);
void backwardSubstitution(int dim, float** U, float* y, float* X);
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


        // Decomposição LU
        clock_t start = clock();
        luDecomposition(DIMENSAO, A, B);
        clock_t end = clock();
        double luTime = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf("Tempo de execucao da Decomposicao LU: %f segundos\n", luTime);
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

void luDecomposition(int dim, float** A, float* B) {
    float** L = (float**)malloc(dim * sizeof(float*));
    float** U = (float**)malloc(dim * sizeof(float*));
    for (int i = 0; i < dim; i++) {
        L[i] = (float*)calloc(dim, sizeof(float));
        U[i] = (float*)calloc(dim, sizeof(float));
    }

    for (int i = 0; i < dim; i++) {
        for (int k = i; k < dim; k++) {
            float soma = 0;
            for (int j = 0; j < i; j++) {
                soma += L[i][j] * U[j][k];
            }
            U[i][k] = A[i][k] - soma;
        }

        for (int k = i; k < dim; k++) {
            if (i == k) {
                L[i][i] = 1;
            } else {
                float soma = 0;
                for (int j = 0; j < i; j++) {
                    soma += L[k][j] * U[j][i];
                }
                L[k][i] = (A[k][i] - soma) / U[i][i];
            }
        }
    }

    float* y = (float*)malloc(dim * sizeof(float));
    float* X = (float*)malloc(dim * sizeof(float));
    forwardSubstitution(dim, L, B, y);
    backwardSubstitution(dim, U, y, X);

    free(y);
    for (int i = 0; i < dim; i++) {
        free(L[i]);
        free(U[i]);
    }
    free(L);
    free(U);

    imprimirSolucoes(dim, X, "SUBSTITUICAO LU");
    free(X);
}

void forwardSubstitution(int dim, float** L, float* B, float* y) {
    for (int i = 0; i < dim; i++) {
        y[i] = B[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
    }
}

void backwardSubstitution(int dim, float** U, float* y, float* X) {
    for (int i = dim - 1; i >= 0; i--) {
        X[i] = y[i];
        for (int j = i + 1; j < dim; j++) {
            X[i] -= U[i][j] * X[j];
        }
        X[i] /= U[i][i];
    }
}

void imprimirSolucoes(int dim, float* X, const char* nomeMetodo) {
    printf("\nSolucoes do sistema pelo metodo %s:\n", nomeMetodo);
    for (int i = 0; i < dim; i++) {
        printf("x%d = %.2f\n", i + 1, X[i]);
    }
}