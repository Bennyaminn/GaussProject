#include <iostream>
#include <fstream>
#include <locale.h>
#include <cmath>

using namespace std;
unsigned wczytaj(double** &A, double* &b);
void gauss_bez_wyboru(double** A, double* b, unsigned matrix_size);
void gauss_z_wyborem_w_kolumnach(double** A, double* b, unsigned matrix_size);
void gauss_z_wyborem_w_wierszach(double** A, double* b, unsigned matrix_size);
void gauus_z_wyborem_pelnym(double** A, double* b, unsigned matrix_size);
void odejmowanie_wiersza(double** A, double* b, unsigned matrix_size, unsigned wiersz_odejmowany, unsigned wiersz_odejmujacy, double mnoznik);
void postepowanie_odwrotne(double** A, double* b, unsigned matrix_size, unsigned* tab_przeksztalcen = nullptr);
void zamiana_wierszy(double** A, double* b, unsigned matrix_size, unsigned zamieniany1, unsigned zamieniany2);
void zamiana_kolumn(double** A, unsigned matrix_size, unsigned zamieniany1, unsigned zamieniany2);
void kopiowanie_tablic(double** A, double* b, unsigned matrix_size, double** &kopia_2D, double* &kopia1D);

int main()
{
    setlocale(0,"");
    int wybor;
    double** A = nullptr;
    double* b = nullptr;
    unsigned matrix_size;

    while(true)
    {
        cout<<"----------MENU----------"<<endl;
        cout<<"1.  Wczytanie danych z pliku."<<endl;
        cout<<"2.  Metoda eliminacji Gaussa bez wyboru elementu."<<endl;
        cout<<"31. Metoda eliminacji Gaussa z wyborem elementu w wierszach."<<endl;
        cout<<"32. Metoda eliminacji Gaussa z wyborem elementu w kolumnach."<<endl;
        cout<<"4.  Metoda eliminacji Gaussa z pe³nym wyborem elementu."<<endl;
        cout<<"5.  Wyjœcie z programu."<<endl<<endl;
        cout<<"Wybór: ";
        cin>>wybor;

        switch(wybor){
        case 1:
            if(matrix_size = wczytaj(A,b))
                cout<<"Wczytano dane!"<<endl;
            break;
        case 2:
            if(A!=nullptr)
                gauss_bez_wyboru(A,b,matrix_size);
            else
                cout<<"Nale¿y najpierw wczytaæ dane!\a"<<endl;
            break;
        case 31:
            if(A!=nullptr)
                gauss_z_wyborem_w_wierszach(A,b,matrix_size);
            else
                cout<<"Nale¿y najpierw wczytaæ dane!\a"<<endl;
            break;
        case 32:
            if(A!=nullptr)
                gauss_z_wyborem_w_kolumnach(A,b,matrix_size);
            else
                cout<<"Nale¿y najpierw wczytaæ dane!\a"<<endl;
            break;
        case 4:
            if(A!=nullptr)
                gauus_z_wyborem_pelnym(A,b,matrix_size);
            else
                cout<<"Nale¿y najpierw wczytaæ dane!\a"<<endl;
            break;
        case 5:
            if(A!=nullptr){
                delete [] b;
                delete [] A[0];
                delete [] A;
            }
            exit(0);
            break;
        default:
            cout<<"Niepoprawny wybór!"<<endl;
        }

        system("pause");
        system("cls");
    }
}

unsigned wczytaj(double** &A, double* &b)
{
    unsigned matrix_size;
    ifstream source_file("dane.csv");
    if (!source_file.is_open())
    {
        cout <<"The file has not been open!"<<endl;
        return 0;
    }
    source_file >> matrix_size;

    A = new double*[matrix_size];
    A[0] = new double[matrix_size*matrix_size];
    for(unsigned i = 1; i< matrix_size; i++)
        A[i] = A[i-1] + matrix_size;

    b = new double[matrix_size];

    char semicolumn;
    for (unsigned i = 0; i < matrix_size+1; i++)
        source_file >> semicolumn;

    for (unsigned i = 0; i < matrix_size; i++)
    {
        for (unsigned j = 0; j < matrix_size; j++)
        {
            source_file >> A[i][j];
            source_file >> semicolumn;
        }
        source_file >> semicolumn;
        source_file >> b[i];
    }
    source_file.close();

    return matrix_size;
}

void gauss_bez_wyboru(double** A, double* b, unsigned matrix_size)
{
    double** A_rob;
    double* b_rob;

    kopiowanie_tablic(A,b,matrix_size, A_rob, b_rob);

    for(unsigned j = 0; j<matrix_size; j++){
        for(unsigned i = j+1; i<matrix_size; i++){
            odejmowanie_wiersza(A_rob, b_rob, matrix_size, i, j, A_rob[i][j]/A_rob[j][j]);
        }
    }

    postepowanie_odwrotne(A_rob, b_rob, matrix_size);

    delete [] b_rob;
    delete [] A_rob[0];
    delete [] A_rob;
}

void odejmowanie_wiersza(double** A, double* b, unsigned matrix_size, unsigned wiersz_odejmowany, unsigned wiersz_odejmujacy, double mnoznik)
{
    for(unsigned i = 0; i<matrix_size; i++){
        A[wiersz_odejmowany][i] -= mnoznik*A[wiersz_odejmujacy][i];
        if(abs(A[wiersz_odejmowany][i]) <= 1e-13) A[wiersz_odejmowany][i]=0;
    }
    b[wiersz_odejmowany] -= mnoznik*b[wiersz_odejmujacy];
}

void postepowanie_odwrotne(double** A, double* b, unsigned matrix_size, unsigned* tab_przeksztalcen)
{
    double x[matrix_size];

    for(int i = matrix_size-1; i>=0; i--){
        for(unsigned j = matrix_size-1; j>i; j--){
                b[i] -= A[i][j]*x[j];
        }
        x[i]=b[i]/A[i][i];
        if(abs(x[i]) <= 1e-13) x[i]=0;
    }

    if(tab_przeksztalcen!=nullptr)
    {
        for(unsigned i = 0; i<matrix_size-1; i++)
            if(tab_przeksztalcen[i]!=i)
                for(unsigned j = i+1; j<matrix_size; j++)
                    if(tab_przeksztalcen[j]==i){
                        swap(tab_przeksztalcen[j],tab_przeksztalcen[i]);
                        swap(x[i],x[j]);
                    }
    }

    cout<<"---WYNIK---"<<endl;
    for(unsigned i = 0; i<matrix_size; i++)
        cout<<"x"<<i+1<<"="<<x[i]<<endl;
    cout<<endl;
}

void gauss_z_wyborem_w_kolumnach(double** A, double* b, unsigned matrix_size)
{
    double** A_rob;
    double* b_rob;

    kopiowanie_tablic(A,b,matrix_size, A_rob, b_rob);

    unsigned max_indeks;
    for(unsigned j = 0; j<matrix_size; j++){
        max_indeks=j;
        for(unsigned k = j+1; k<matrix_size; k++)
            if(A_rob[k][j]>A_rob[max_indeks][j])
                max_indeks=k;

        if(max_indeks!=j)
            zamiana_wierszy(A_rob,b_rob,matrix_size,max_indeks,j);

        for(unsigned i = j+1; i<matrix_size; i++)
            odejmowanie_wiersza(A_rob, b_rob, matrix_size, i, j, A_rob[i][j]/A_rob[j][j]);
    }

    postepowanie_odwrotne(A_rob, b_rob, matrix_size);

    delete [] b_rob;
    delete [] A_rob[0];
    delete [] A_rob;
}

void zamiana_wierszy(double** A, double* b, unsigned matrix_size, unsigned zamieniany1, unsigned zamieniany2)
{
    for(unsigned i = 0; i<matrix_size; i++)
    {
        swap(A[zamieniany1][i],A[zamieniany2][i]);
    }
    swap(b[zamieniany1],b[zamieniany2]);
}

void kopiowanie_tablic(double** A, double* b, unsigned matrix_size, double** &kopia_2D, double* &kopia_1D)
{
    kopia_2D = new double*[matrix_size];
    kopia_2D[0] = new double[matrix_size*matrix_size];
    for(unsigned i = 1; i< matrix_size; i++)
        kopia_2D[i] = kopia_2D[i-1] + matrix_size;
    kopia_1D = new double[matrix_size];

    for(int i=0;i<matrix_size;i++){
        for(int j=0;j<matrix_size;j++){
            kopia_2D[i][j]=A[i][j];
        }
        kopia_1D[i]=b[i];
    }
}

void gauss_z_wyborem_w_wierszach(double** A, double* b, unsigned matrix_size)
{
    double** A_rob;
    double* b_rob;

    kopiowanie_tablic(A,b,matrix_size, A_rob, b_rob);

    unsigned max_indeks;
    unsigned tablica_przeksztalcen[matrix_size];
    for(unsigned i = 0; i<matrix_size; i++)
        tablica_przeksztalcen[i]=i;

    for(unsigned j = 0; j<matrix_size; j++)
    {
        max_indeks=j;
        for(unsigned k = j+1; k<matrix_size; k++)
            if(A_rob[j][k]>A_rob[j][max_indeks])
                max_indeks=k;
        if(max_indeks!=j){
                zamiana_kolumn(A_rob,matrix_size,j,max_indeks);
                swap(tablica_przeksztalcen[j],tablica_przeksztalcen[max_indeks]);
        }

        for(unsigned i = j+1; i<matrix_size; i++)
            odejmowanie_wiersza(A_rob, b_rob, matrix_size, i, j, A_rob[i][j]/A_rob[j][j]);
    }

    postepowanie_odwrotne(A_rob, b_rob, matrix_size, &tablica_przeksztalcen[0]);

    delete [] b_rob;
    delete [] A_rob[0];
    delete [] A_rob;

}

void zamiana_kolumn(double** A, unsigned matrix_size, unsigned zamieniany1, unsigned zamieniany2)
{
    for(unsigned i = 0; i<matrix_size; i++)
    {
        swap(A[i][zamieniany1],A[i][zamieniany2]);
    }
}

void gauus_z_wyborem_pelnym(double** A, double* b, unsigned matrix_size)
{
    double** A_rob;
    double* b_rob;

    kopiowanie_tablic(A,b,matrix_size, A_rob, b_rob);

    unsigned max_indeks_wiersza, max_indeks_kolumny;
    unsigned tablica_przeksztalcen[matrix_size];
    for(unsigned i = 0; i<matrix_size; i++)
        tablica_przeksztalcen[i]=i;

    for(unsigned j = 0; j<matrix_size; j++)
    {
        max_indeks_wiersza=j;
        max_indeks_kolumny=j;

        for(unsigned w = j; w<matrix_size; w++)
            for(unsigned k = j; k<matrix_size; k++)
            {
                if(A_rob[w][k]>A_rob[max_indeks_wiersza][max_indeks_kolumny])
                {
                    max_indeks_wiersza=w;
                    max_indeks_kolumny=k;
                }
            }

        if(max_indeks_kolumny!=j){
                zamiana_kolumn(A_rob,matrix_size,j,max_indeks_kolumny);
                swap(tablica_przeksztalcen[j],tablica_przeksztalcen[max_indeks_kolumny]);
        }
        if(max_indeks_wiersza!=j)
            zamiana_wierszy(A_rob,b_rob,matrix_size,max_indeks_wiersza,j);

        for(unsigned i = j+1; i<matrix_size; i++)
            odejmowanie_wiersza(A_rob, b_rob, matrix_size, i, j, A_rob[i][j]/A_rob[j][j]);
    }

    postepowanie_odwrotne(A_rob, b_rob, matrix_size, &tablica_przeksztalcen[0]);

    delete [] b_rob;
    delete [] A_rob[0];
    delete [] A_rob;
}
