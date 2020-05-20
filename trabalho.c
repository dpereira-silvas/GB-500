#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


void E( mpz_t c, mpz_t m, mpz_t g, mpz_t n ) {
    mpz_t r;
    mpz_t n2;
    
    mpz_init( r );
    mpz_init( n2 );
    
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, time(NULL) );
    
    mpz_mul( n2, n, n );
    
    mpz_urandomm( r, rstate, n );
    mpz_powm( c, g, m, n2 );
    mpz_powm( r, r, n, n2 );
    mpz_mul( c, c, r );    
    mpz_mod( c, c, n2 );     // c1 = (g^m1)(r^n) mod n^2
}

void D( mpz_t m, mpz_t c, mpz_t lambda, mpz_t micro, mpz_t n ) {
    mpz_t r;
    mpz_t n2;
    mpz_t L;
    
    mpz_init( n2 );
    mpz_init( L );
    
    mpz_mul( n2, n, n );
    
    mpz_powm( L, c, lambda, n2 );
    mpz_sub_ui( L, L, 1 );
    mpz_div( L, L, n );
    
    mpz_mul( m, L, micro );
    mpz_mod( m, m, n );
    
}

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
        mpz_set(u[i], v[i]);
    }
}


int main() {
    int nbits = 2048;
    
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, time(NULL) );
    
    
    mpz_t p;
    mpz_t q;
    mpz_t n;
    mpz_t n2;
    mpz_t g;
    mpz_t lambda;
    mpz_t micro;
    
    
    mpz_init( p );
    mpz_init( q );
    mpz_init( n );
    mpz_init( n2 );
    mpz_init( g );
    mpz_init( lambda );
    mpz_init( micro );
    
    mpz_urandomb( p, rstate, nbits/2 );
    mpz_urandomb( q, rstate, nbits/2 );
    
    mpz_nextprime( p, p );
    mpz_nextprime( q, q );
        
    mpz_mul( n, p, q );    // n = p * q

    mpz_sub_ui( p, p, 1 ); // p = p - 1
    mpz_sub_ui( q, q, 1 ); // q = q - 1
    
    mpz_mul( lambda, p, q ); // lambda = (p-1)(q-1)
    
    
    mpz_invert( micro, lambda, n );  // micro = lambda^-1 mod n
    
    mpz_add_ui( g, n , 1 );  // g = n + 1 
    
    mpz_mul( n2, n, n );
    
    // publica : (n,g)
    // privada : (lambda, micro)
    

    FILE *arq;
    arq = fopen("out.txt","w");
    
    int dx = 5;
    int dt = 5;
    // Exemplo 1
    // int n_nodes = 11;

    // Exemplo 2
    int n_nodes = 7;

    int t_final = 500;

    // numero de Fourier 1/mul
    int mul = 5;

    mpz_t *Uant;
    Uant = malloc(sizeof(mpz_t) * n_nodes);

    // Condições de contorno
    // Exemplo 1
    // mpz_init_set_str(Uant[0],"0",10);
    // mpz_init_set_str(Uant[n_nodes-1],"0",10);

    // Exemplo 2
    mpz_init_set_str(Uant[0],"20000000",10);
    mpz_init_set_str(Uant[n_nodes-1],"50000000",10);

    // Condição inicial
    for(int i = 1; i < (n_nodes-1); i++)
    {   
        //Exemplo 1
        // mpz_init_set_str(Uant[i],"20000000",10); // 20.000000

        // Exemplo 2
        mpz_init_set_ui(Uant[i],(60-2*(i*dx))*1000000  ); // 60 -2*x
    }
    // print(Uant,n_nodes);
    write_vector(Uant, n_nodes-1, arq);


    mpz_t *U;
    U = malloc(sizeof(mpz_t) * n_nodes);

    // inicializando vetor
    for(int i = 0; i < n_nodes; i++)
    {
        mpz_init(U[i]); 
    }
    


    mpz_t *C;
    C = malloc(sizeof(mpz_t) * n_nodes);
    // inicializando vetor
    for(int i = 0; i < n_nodes; i++)
    {
        mpz_init(C[i]); 
    }

    mpz_t *U_dec;
    U_dec = malloc(sizeof(mpz_t) * n_nodes);
    for(int i = 0; i < n_nodes; i++)
    {
        mpz_init(U_dec[i]); 
    }



    // Cifrando a condição inicial
    for(int i = 0; i < n_nodes; i++)
    {
        E( C[i], Uant[i], g, n );
    }
    clear_vector(Uant,n_nodes);

    // passando para o vetor U as condições de contorno criptografadas
    mpz_set(U[0],C[0]);
    mpz_set(U[n_nodes-1],C[n_nodes-1]);

    mpz_t tmp;
    mpz_init(tmp);

    int j = dt; // Iteração começa do instante de tempo 0+dt

    while(j <= t_final)
    {
        /////   SERVIDOR
        for(int i = 1; i < (n_nodes-1); i++)
        {
            mpz_powm_ui( tmp, C[i], 2, n2 );      // 
            mpz_invert(tmp, tmp, n2);            //  -2*Uant[i]

            mpz_powm_ui( U[i], C[i], mul, n2 );  // 5*Uant[i]

            mpz_mul( U[i], U[i], C[i+1] );       //  U[i] = U[i] + Uant[i+1]
            mpz_mod( U[i], U[i], n2 );

            mpz_mul( U[i], U[i], C[i-1] );       //  U[i] = U[i] + Uant[i-1]
            mpz_mod( U[i], U[i], n2 );

            mpz_mul( U[i], U[i], tmp );          //  U[i] = U[i] -  2*Uant[i]
            mpz_mod( U[i], U[i], n2 );
        }


        /////   CLIENTE
        // Decifra o U
        for(int i = 0; i < n_nodes; i++)
        {
            D( U_dec[i], U[i], lambda, micro, n );
            if((i != 0) && (i != (n_nodes-1)) )
                mpz_div_ui(U_dec[i],U_dec[i],mul); 
        }

        // print(U_dec,n_nodes);
        write_vector(U_dec, n_nodes, arq);

        // Cifra novamente e manda para o servidor
        for(int i = 0; i < n_nodes; i++)
        {
            E( C[i], U_dec[i], g, n );
        }
        j = j + dt;
    }

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(n);
    mpz_clear(n2);
    mpz_clear(g);
    mpz_clear(lambda);
    mpz_clear(micro);

    mpz_clear(tmp);
    clear_vector(U,n_nodes);
    clear_vector(C,n_nodes);
    clear_vector(U_dec,n_nodes);
    fclose(arq);

    return 0;
}
