#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void print(mpz_t *v, int n)
{
    for(int i = 0; i < (n-1); i++)
    {
        printf("%f ",(1/1000000.0)*mpz_get_ui(v[i] ));
    }
    printf("%f\n",(1/1000000.0)*mpz_get_ui(v[(n-1)]));
}

void write_vector(mpz_t *v, int n, FILE *f)
{

    for(int i = 0; i < (n-1); i++)
    {
        fprintf (f, "%f\t",(1/1000000.0)*mpz_get_ui(v[i] ));
    }
    fprintf(f,"%f\n",(1/1000000.0)*mpz_get_ui(v[(n-1)]));  
}



void clear_vector(mpz_t *v, int n)
{
    for(int i = 0; i < n; i++)
    {
        mpz_clear(v[i]);
    }
        free(v);

}

void copy_vector(mpz_t *u, mpz_t *v, int n)
{
    for(int i = 0; i < n; i++)
    {
        // u[i] = v[i];
        mpz_set(u[i], v[i]);
    }

}


int main() {
        

    FILE *arq;
    arq = fopen("NH_out.txt","w");
    
    int dx = 5;
    int dt = 5;
    // int n_nodes = 11;
    int n_nodes = 7;
    int t_final = 500;
    


    // numero de Fourier 1/mul
    int mul = 5;


    mpz_t *Uant;
    Uant = malloc(sizeof(mpz_t) * n_nodes);
    // Condições de contorno
    // mpz_init_set_str(Uant[0],"0",10);
    // mpz_init_set_str(Uant[n_nodes-1],"0",10);

    mpz_init_set_str(Uant[0],"20000000",10);
    mpz_init_set_str(Uant[n_nodes-1],"50000000",10);

    // Condição inicial
    for(int i = 1; i < (n_nodes-1); i++)
    {
        // mpz_init_set_str(Uant[i],"20000000",10); // 20.000000

        mpz_init_set_ui(Uant[i],(60-2*(i*dx))*1000000  ); // 60 -2*x
    }
    // mpz_init_set_str(Uant[n_nodes-1],"0",10);

    print(Uant,n_nodes);
    write_vector(Uant, n_nodes, arq);

    mpz_t *U;
    U = malloc(sizeof(mpz_t) * n_nodes);

    // inicializando vetor
    for(int i = 0; i < n_nodes; i++)
    {
        mpz_init(U[i]); 
    }
    mpz_set(U[0],Uant[0]);
    mpz_set(U[n_nodes-1],Uant[n_nodes-1]);

    
    
    mpz_t tmp;
    mpz_init(tmp);

    int j = dt; // Iteração começa do instante de tempo 0+dt
    while(j <= t_final)
    {
        /////   SERVIDOR
        for(int i = 1; i < (n_nodes-1); i++)
        {
            
            mpz_mul_ui(tmp,Uant[i],2);          // tmp = 2*Uant[i]
            mpz_div_ui(tmp,tmp,mul);            // tmp = (0.2)*2*Uant[i]

            mpz_add(U[i],Uant[i+1],Uant[i-1]);  // U[i] = Uant[i+1]+Uant[i-1]
            mpz_div_ui(U[i],U[i],mul);          // U[i] = (0.2)*U[i]
            mpz_add(U[i],U[i],Uant[i]);         // U[i] = (0.2)*(Uant[i+1]+Uant[i-1]) + Uant[i]

            mpz_sub(U[i], U[i], tmp);           // U[i] = (0.2)*(Uant[i+1]+Uant[i-1]) + Uant[i] - (0.2)*2*Uant[i]
        }
        print(U,n_nodes);
        copy_vector(Uant,U,n_nodes);
        
        write_vector(U, n_nodes, arq);

        j = j + dt;
    }

    clear_vector(U,n_nodes);
    fclose(arq);

    return 0;
}
