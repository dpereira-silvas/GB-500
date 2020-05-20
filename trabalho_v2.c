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



void print(mpz_t *v, int n, int op)
{
    if(op)
    {
        for(int i = 0; i < (n-1); i++)
        {
            printf("%f ",(1/5000000.0)*mpz_get_ui(v[i] ));
        }
        printf("%f\n",(1/5000000.0)*mpz_get_ui(v[(n-1)]));
    }
    else
    {
            for(int i = 0; i < (n-1); i++)
            {
                printf("%f ",(1/1000000.0)*mpz_get_ui(v[i] ));
            }
            printf("%f\n",(1/1000000.0)*mpz_get_ui(v[(n-1)]));

    }

}

void write_vector(mpz_t *v, int n, FILE *f, int op)
{
    if(op)
    {
        for(int i = 0; i < (n-1); i++)
        {
            fprintf (f, "%f\t",(1/5000000.0)*mpz_get_ui(v[i] ));
        }
        fprintf(f,"%f\n",(1/5000000.0)*mpz_get_ui(v[(n-1)]));
    }
    else
    {
        for(int i = 0; i < (n-1); i++)
        {
            fprintf (f, "%f\t",(1/1000000.0)*mpz_get_ui(v[i] ));
        }
        fprintf(f,"%f\n",(1/1000000.0)*mpz_get_ui(v[(n-1)]));
    }
    
    
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
    int t_final = 500;
    float l = 0.2;

    mpz_t *Uant;
    Uant = malloc(sizeof(mpz_t) * 11);
    mpz_init_set_str(Uant[0],"0",10);
    mpz_init_set_str(Uant[10],"0",10);

    // setando a condição inicial
    for(int i = 1; i < 10; i++)
    {
        mpz_init_set_str(Uant[i],"20000000",10); // 20.000000
    }
    print(Uant,11,0);
    write_vector(Uant, 11, arq,0);


    mpz_t *U;
    U = malloc(sizeof(mpz_t) * 11);
    // inicializando vetor
    for(int i = 0; i < 11; i++)
    {
        mpz_init(U[i]); 
    }

    // print(Uant,11);
    // write_vector(Uant, 11, arq);
    
    int j = dt;
    
    mpz_t *C;
    C = malloc(sizeof(mpz_t) * 11);
    // inicializando vetor
    for(int i = 0; i < 11; i++)
    {
        mpz_init(C[i]); 
    }

    mpz_t *U_dec;
    U_dec = malloc(sizeof(mpz_t) * 11);
    for(int i = 0; i < 11; i++)
    {
        mpz_init(U_dec[i]); 
    }

    // Cifrando a condição inicial
    for(int i = 0; i < 11; i++)
    {
        E( C[i], Uant[i], g, n );
    }
    clear_vector(Uant,11);

    mpz_set(U[0],C[0]);
    mpz_set(U[10],C[10]);

    mpz_t mul;
    mpz_init_set_str( mul, "5", 10 );
    mpz_invert(mul, mul, n2);


    mpz_t tmp1;
    mpz_init(tmp1);
    // mpz_t tmp2;
    // mpz_init(tmp2);


    while(j <= t_final)
    {
        /////   SERVIDOR
        for(int i = 1; i < 10; i++)
        {
            mpz_powm_ui( tmp1, C[i], 2, n2 );  // 
            mpz_invert(tmp1, tmp1, n2);        //  -2*Uant[i] mod n2
            mpz_powm( tmp1, tmp1, mul, n2 );  // -(2/5)*Uant[i] mod n2

            // mpz_powm_ui( tmp2, C[i], 5, n2 );  // 5*Uant[i]

            mpz_set(U[i],C[i+1]);              //  U[i] = Uant[i+1]
            mpz_mul( U[i], U[i], C[i-1] );     //  U[i] = U[i] + Uant[i-1] mod n2
            mpz_mod( U[i], U[i], n2 );
            mpz_powm( U[i], U[i], mul, n2 );  // -(1/5)*U[i] mod n2
            
            mpz_mul( U[i], U[i], C[i] );       //  U[i] = U[i] + Uant[i] mod n2
            mpz_mod( U[i], U[i], n2 );
            
            
            if(mpz_cmp(U[i],tmp1) < 0)
            {
                printf("Deu merda %d\n",i);

            }
            mpz_mul( U[i], U[i], tmp1 );       //  U[i] = U[i] -  (2/5)*Uant[i] mod n2
            mpz_mod( U[i], U[i], n2 );

       
            // mpz_powm( U[i], U[i], mul, n2 );    // U[i] = (1/mul)*U[i] mod n2

            // mpz_mul( U[i], U[i], C[i] );       //  U[i] = U[i] + Uant[i] mod n2
            // mpz_mod( U[i], U[i], n2 );

        }


        /////   CLIENTE
        // decifra o U
        for(int i = 0; i < 11; i++)
        {
            D( U_dec[i], U[i], lambda, micro, n );
            // mpz_div_ui(U_dec[i],U_dec[i],5); 
        }

        print(U_dec,11,0);
        write_vector(U_dec, 11, arq,0);
        // break;
        copy_vector(C, U, 11);
        // Divide todo o U por 5
        // for(int i = 1; i < 10; i++)
        // {
        //     mpz_div_ui(U_dec[i],U_dec[i],5); 
        // }

        // cifra novamente e manda para o servidor
        // for(int i = 0; i < 11; i++)
        // {
        //     E( C[i], U_dec[i], g, n );
        // }
        j = j + dt;
    }

    clear_vector(U,11);
    clear_vector(C,11);
    clear_vector(U_dec,11);
    fclose(arq);

    return 0;
}
