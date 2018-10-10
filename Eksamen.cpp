#define NMAX 20

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>

using namespace std;

//---------------------------------------------------------------------------------------
//------------------------- DEKLARATION AF ANVENDTE FUNKTIONER --------------------------
//---------------------------------------------------------------------------------------

//----------------- FUNKTIONER TIL INDHENTNING OG UDSKRIVNING AF DATA -------------------
//Indhentning:
void AngivM(double M[NMAX][NMAX], int n, int m);
void Angivv(double v[NMAX], int n);
void HentM(double M[NMAX][NMAX], const char* filename);

//Udskrivning:
void UdskrivM(double M[NMAX][NMAX], int n, int m);
void Udskrivv(double v[NMAX], int n);
void GemBogC(double A[NMAX][NMAX], double b[NMAX], int n, int m, const char* filename);

//----------------- FUNKTIONER TIL INTRODUKTIONEN ---------------------------------------
void Introduktion();
void CaseA(double A[NMAX][NMAX], double b[NMAX], int &n, int &m);
void CaseB(double A[NMAX][NMAX], double b[NMAX], int &n, int &m);
void BeregnAogb(double A[NMAX][NMAX], double b[NMAX], int &n, int &m);
void CaseC(double A[NMAX][NMAX], double b[NMAX], int &n, int &m);

//----------------- FUNKTIONER TIL MATRIX OG VEKTOR REGNING -----------------------------
//Matrix:
void ProduktMogv(double A[NMAX][NMAX], double x[NMAX], double b[NMAX], int n, int m);
void ProduktMogM(double A[NMAX][NMAX],double B[NMAX][NMAX], double P[NMAX][NMAX], int n, int m);
void TransponerM(double A[NMAX][NMAX],double AT[NMAX][NMAX],int n, int m);
void ProduktMogMT(double A[NMAX][NMAX],double AT[NMAX][NMAX], double P[NMAX][NMAX], int n, int m);
void Totalmatrix(double TM[NMAX][NMAX+1],double A[NMAX][NMAX],double y[NMAX],int m);
void Tagsoejle(double A[NMAX][NMAX], double v[NMAX], int n, int k);
void RetMatrix(double M[NMAX][NMAX],int n, int m);

//Vektor:
void Indprodukt(double a[NMAX], double b[NMAX], double p[NMAX], int n);
double Pprodukt(double a[NMAX], double b[NMAX], int n);
double Ldiffvek(double a[NMAX], double b[NMAX], double n);
double Lvek(double v[NMAX], int n);
void RetVektor(double v[NMAX],int n);

//----------------- FUNKTIONER TIL DE ENKELTE METODER -----------------------------------
//Vigtige funktioner:
void BackwardsSubstitution(double TM[NMAX][NMAX+1],double x[NMAX],int n);

//Metode 1:
void DelvisPivotering(double TM[NMAX][NMAX+1],int j,int n);
void Gauss(double TM[NMAX][NMAX+1],int n,int &bs);
void Metode1(double A[NMAX][NMAX],double b[NMAX],double x[NMAX],int n, int m);

//Metode 2:
void DanQ(double A[NMAX][NMAX], double Q[NMAX][NMAX], double R[NMAX][NMAX], int n, int k);
void QRmodGS(double A[NMAX][NMAX], double Q[NMAX][NMAX], double R[NMAX][NMAX], int n, int m);
void Metode2(double A[NMAX][NMAX], double b[NMAX], double Q[NMAX][NMAX], double R[NMAX][NMAX], int n, int m);

//Metode 3:
void ParameterMetode3(double x[NMAX], double &eps,  int &N, int m);
double f(double A[NMAX][NMAX], double b[NMAX], double x[NMAX], int n, int m);
double falpha(double A[NMAX][NMAX], double b[NMAX], double x[NMAX], double sk[NMAX], double alpha, int n, int m);
void Gradient(double A[NMAX][NMAX], double b[NMAX], double x[NMAX], double grad[NMAX], int n, int m);
void DanHesse(double H[NMAX][NMAX], int m);
void Hesse(double A[NMAX][NMAX], double b[NMAX], double H[NMAX][NMAX], double xny[NMAX], double xgl[NMAX], int n, int m);
void BFGS(double H[NMAX][NMAX], double grad[NMAX], double sk[NMAX], int m);
void IndkredsningsAlg(double A[NMAX][NMAX], double vb[NMAX], double x[NMAX], double sk[NMAX], double d, int n, int m, double &a, double &b);
void GoldenSection(double A[NMAX][NMAX], double b[NMAX], double x[NMAX],double sk[NMAX],int n,int m,double &alph);
void Metode3(double A[NMAX][NMAX], double H[NMAX][NMAX], double b[NMAX], int n, int m);

//Analyse af Metode 3:
void Danenhedsvek(double e[NMAX], int m, int k);
void DanRinv(double R[NMAX][NMAX], double Rinv[NMAX][NMAX], int m);
void RinvogRinvT(double Rinv[NMAX][NMAX], double P[NMAX][NMAX], int m);
void AnalyseMetode3(double R[NMAX][NMAX], double P[NMAX][NMAX], double H[NMAX][NMAX], int m);

//---------------------------------------------------------------------------------------
//---------------------------------------- MAIN -----------------------------------------
//---------------------------------------------------------------------------------------

int main() {

    char Valgabc, Valg1eller23, Valg2eller3, ValgIgen, ValgGemB, ValgGemC;
    double A[NMAX][NMAX], Q[NMAX][NMAX], R[NMAX][NMAX];
    double H[NMAX][NMAX], P[NMAX][NMAX];
    double b[NMAX], x[NMAX];
    int n,m,retry;

    do{
        Introduktion();
        do{
            cout<<"\nVælg mellem a, b eller c: ";   cin>>Valgabc;
            retry=0;

            switch (Valgabc)
            {

            case 'a' : {
                    cout<<"\nDu valgte a.";
                    cout<<"\nPå dette sted indlæses n, m, A og B fra fil."<<endl;

                    CaseA(A,b,n,m);
                }
                break;

            case 'b' : {
                    cout<<"\nDu valgte b.";
                    cout<<"\nPå dette sted indtastes n, m, A og B manuelt."<<endl;

                    CaseB(A,b,n,m);

                    cout<<"\nØnsker du at skrive n, m, A og b til fil? J/N: ";    cin>>ValgGemB;
                    if(ValgGemB=='J')
                    {
                        GemBogC(A,b,n,m,"GemB.txt");
                    }
                }
                break;

            case 'c' : {
                    cout<<"\nDu valgte c.";
                    cout<<"\nPå dette sted indtastes n, m og R manuelt."<<endl;
                    cout<<"\nFrembringer matricen Q og beregner A = Q * R.";
                    cout<<"\nIndtastning x* og t*, derefter beregning bp, e og b"<<endl;

                    CaseC(A,b,n,m);

                    cout<<"\nØnsker du at skrive n, m, A og b til fil? J/N: ";  cin>>ValgGemC;
                    if(ValgGemC=='J')
                    {
                        GemBogC(A,b,n,m,"GemC.txt");
                    }
                }
                break;

            default :{
                    cout<<"Du kan kun vælge a, b eller c"<<endl;
                    retry=1;
                }

            }// Switch [SLUT]
        }while(retry==1);
        do{
            retry=0;

            cout<<"\nDu kan nu vælge følgende fremgangsmåder:"<<endl;
            cout<<"\n1) Bestemmelse af x* ved normalligningerne: At * A * x = At * b"<<endl;
            cout<<"\n   Eller efter udførelse af QR-faktorisering af matrix A med modificeret Gram-Schmidt:"<<endl;
            cout<<"\n2) Bestemmelse af x* ved at løse R * x = Qt * b:";
            cout<<"\n3) Bestemmelse af x* ved minimering af f(x) = (b - A * x)^2"<<endl;
            cout<<"\nVælg mellem (1) eller (2/3): ";    cin>>Valg1eller23;

            if(Valg1eller23=='1')
            {
                  Metode1(A, b, x, n, m);
            }

            else if (Valg1eller23 == '2' || Valg1eller23 == '3')
            {
                QRmodGS(A,Q,R,n,m);
                do{
                    cout<<"\nDu har valgt 2 eller 3"<<endl;
                    cout<<"\nQR-faktorisering af matrix A er udført ved hjælp af modificeret Gram-Schmidt til:"<<endl;cout<<"\nMatrix A: "<<endl;
                    UdskrivM(A,n,m);
                    cout<<endl<<"Matrix Q:"<<endl;
                    UdskrivM(Q,n,m);
                    cout<<endl<<"Matrix R:"<<endl;
                    UdskrivM(R,m,m);
                    cout<<"\nDu kan nu vælge imellem:"<<endl;
                    cout<<"\n2) Bestemmelse af mindre kvadraters løsningen x* ved at løse R * x = Qt * b";
                    cout<<"\n3) Bestemmelse af mindre kvadraters løsningen x* ved minimering af f(x) = (b - A * x)^2"<<endl;
                    cout<<"\nVælg mellem 2 eller 3: ";  cin>>Valg2eller3;
                    retry=0;
                    switch (Valg2eller3) // Switch [START]
                    {
                        case '2':
                        {
                            cout<<"\nDu valgte 2, der løses efter ligningen R * x = Qt * b"<<endl;
                            Metode2(A,b,Q,R,n,m); //Udskriver QR-faktorisering
                        }
                            break;

                        case '3':
                        {
                            cout<<"\nDu valgte 3, der løses efter minimering af f(x) = (b - Ax)^2"<<endl;
                            Metode3(A,H,b,n,m);
                            AnalyseMetode3(R,P,H,m);
                        }
                            break;
                        default :
                        {
                            cout<<"nej, du skal skriv 2 eller 3"<<endl;
                            retry=1;
                        }
                    } //Switch [SLUT]
                }while(retry==1);
            }
            else{
            	cout<<"\nDu kan kun skrive 1, 2 eller 3, prøv igen: "<<endl;
            	retry=1;
            }
        }while(retry == 1);

        cout << "\nSkal programmet startes fra begyndelsen igen? Tast J/N: ";
        cin >> ValgIgen;

    }while(ValgIgen=='J');
}

// --------------------------------------------------------------------------------------
// --------------------------------- KODE TIL FUNKTIONER --------------------------------
// --------------------------------------------------------------------------------------

//----------------- FUNKTIONER TIL INDHENTNING OG UDSKRIVNING AF DATA -------------------

//Indhentning:
void AngivM(double M[NMAX][NMAX], int n, int m)
{
    for(int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            cout << "Angiv elementet [" << i+1 << "," << j+1 << "] = ";
            cin>>M[i][j];
        }
    }
}
void Angivv(double v[NMAX], int n)
{
    for(int i=0;i<n;i++){
        cout << "Angiv " << i+1 << ". element: ";
        cin>>v[i];
    }
}
void HentM(double M[NMAX][NMAX], const char* filename)
{
    ifstream Fil;
    Fil.open(filename);
    for (int i=0;i<8;i++){
        for (int j=0;j<8;j++) Fil >> M[i][j];
    }
    Fil.close();
}

//Udskrivning:
void UdskrivM(double M[NMAX][NMAX], int n, int m)
{
    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++) cout << setw(10) << M[i][j];
        cout << endl;
    }
}
void Udskrivv(double v[NMAX], int n)
{
    cout << setw(4) << "[" << v[0];
    for (int i=1;i<n;i++) cout << ", " <<v[i];
    cout << "]" << endl;
}
void GemBogC(double A[NMAX][NMAX], double b[NMAX], int n, int m, const char* filename)
{
    ofstream Fil;
    Fil.open(filename);
    Fil << n << endl;
    Fil << m << endl;
    for (int i=0;i<n;i++){
        for (int j=0;j<=m;j++){
            if(j<m) Fil << A[i][j] << " ";
            else Fil << b[i] << endl;
        }
    }
    Fil.close();
}

//----------------- FUNKTIONER TIL INTRODUKTIONEN ---------------------------------------
void Introduktion()
{
    cout<<"Velkommen til programmet!"<<endl<<endl;
    cout<<"Formål: Definition af testproblem A * x = b med kendt"<<endl;
    cout<<"        mindste kvadraters løsning x* og e = b - A * x* "<<endl<<endl;
    cout<<"Du har følgende muligheder for tilvejebringelse af n, m, A og b:"<<endl;
    cout<<"\na) Indlæsning af data fra fil.";
    cout<<"\nb) Indtastning af n, m, A og b, med mulighed for udskrift til fil.";
    cout<<"\nc) Dannelse af Q og indtastning af R til beregning af A = Q * R";
    cout<<"\n   Beregning af af bp og e ud fra indtastet x* og t* samt beregning af b = bp + e"<<endl;
}
void CaseA(double A[NMAX][NMAX],double b[NMAX], int &n, int &m){
    int i,j;

    ifstream IndFil;

    IndFil.open("HentFil.txt");

    IndFil>>n;
    IndFil>>m;

	for (i=0;i<n;i++){
		for (j=0;j<=m;j++){
			if(j<m) IndFil >> A[i][j];
			else IndFil >> b[i];
		}
	}
	IndFil.close();

    cout<<"\nn = "<<n;
    cout<<"\nm = "<<m<<endl;

    cout<<"\nMatrix A: "<<endl;
    UdskrivM(A,n,m);

    cout<<"\nVektor b: "<<endl;
    Udskrivv(b,n);

    RetMatrix(A,n,m);
    RetVektor(b,n);
}
void CaseB(double A[NMAX][NMAX], double b[NMAX], int &n, int &m)
{
    cout<<"\nMatricen er en (n x m) matrix hvor n > m"<<endl;
    cout<<"\nIndtast n: ";  cin>>n;
    cout<<"\nIndtast m: ";  cin>>m;
    AngivM(A,n,m);
    cout<<"\nMatrix A er givet ved:"<<endl;
    UdskrivM(A,n,m);
    RetMatrix(A,n,m);

    cout<<"\nIndtast vektor b: "<<endl;
    Angivv(b,n);
    cout<<"\nVektor b er givet ved:"<<endl;
    Udskrivv(b,n);
    RetVektor(b,n);
}
void BeregnAogb(double A[NMAX][NMAX], double b[NMAX], int &n, int &m)
{
    double D[NMAX][NMAX], M[NMAX][NMAX],Q[NMAX][NMAX],C[NMAX][NMAX], R[NMAX][NMAX];
    double x[NMAX],v[NMAX], bp[NMAX], t[NMAX], e[NMAX];

    HentM(M,"DelC.txt");
    cout << "\nAngiv om du vil have (4) eller (8) rækker, n = ";
    cin >> n;
    cout << "\nMatrix A"<<n<<" er givet ved: "<<endl;
    UdskrivM(M,n,n);

	for(int i=0;i<n;i++){
		Tagsoejle(M,v,n,i);
		D[i][i]=1/sqrt(Pprodukt(v,v,n));
	}
    cout << "\nMatrix D"<<n<<" er givet ved: "<<endl;
	UdskrivM(D,n,n);
    cout << "\nAngiv antallet af søjler hvor m < n, m = ";
    cin >> m;
	ProduktMogMT(D,M,Q,n,n);
    cout << "\nMatrix Q er givet ved: "<<endl;
    UdskrivM(Q,n,m);
    cout << "\nAngiv R's elementer:"  << endl;
    for (int i=0;i<m;i++){
        for (int j=0;j<m;j++){
            if (j>=i){
                cout << "Angiv elementet [" << i+1 << "," << j+1 << "] = ";
                cin >> R[i][j];
            }
            else R[i][j] = 0;
        }
    }
    ProduktMogM(Q,R,A,n,m);
    cout << "\nAngiv elementerne i x*" << endl;
    Angivv(x,m);
    ProduktMogv(A,x,bp,n,m);
    cout << "\nAngiv elementerne i t*" << endl;
    Angivv(t,(n-m));
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if (j>=m) C[i][(j-m)]=M[i][j];
        }
    }
    ProduktMogv(C,t,e,n,(n-m));
    for (int i=0;i<n;i++)
    {
        b[i] = bp[i] + e[i];
    }
    cout << "\nOrtogonalprojektionen bp er givet ved:"<<endl;
    Udskrivv(bp,n);
    cout << "\nFejlvektoren e er givet ved:"<<endl;
    Udskrivv(e,n);
    cout << endl;
}
void CaseC(double A[NMAX][NMAX], double b[NMAX], int &n, int &m)
{
    BeregnAogb(A,b,n,m);
    cout << "Matrix A er givet ved:" << endl;
    UdskrivM(A,n,m);
    cout << "\nVektor b er givet ved:"<<endl;
    Udskrivv(b,n);

    RetMatrix(A,n,m);
    RetVektor(b,n);
}

//----------------- FUNKTIONER TIL MATRIX OG VEKTOR REGNING -----------------------------
//Matrix:
void RetMatrix(double M[NMAX][NMAX],int n, int m)
{
    int i,j;
    char Svar='n';

    cout<<"\nSkal der rettes i matricen? J/N: ";
    cin>>Svar;

    while(Svar=='J') {
        cout<<"\nIndtast række nr: ";
        cin>>i;i--;

        cout<<"\nIndtast søjle nr: ";
        cin>>j;j--;

        cout<<"\nIndtast ny værdi: ";
        cin>>M[i][j];

        cout<<"\nNy matrix givet ved: "<<endl;
        UdskrivM(M,n,m);

        cout<<"\nSkal der rettes i matricen igen? J/N: ";
        cin>>Svar;
    }
}
void ProduktMogv(double A[NMAX][NMAX], double x[NMAX], double b[NMAX], int n, int m)
{
    for (int i=0;i<n;i++)
    {
        b[i] = 0;
        for (int j=0;j<m;j++) b[i] += A[i][j]*x[j];
    }
}
void ProduktMogM(double A[NMAX][NMAX],double B[NMAX][NMAX], double P[NMAX][NMAX], int n, int m)
{
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++)
        {
            P[i][j] = 0;
            for (int k=0;k<m;k++) P[i][j] = P[i][j] + A[i][k]*B[k][j];
        }
    }
}
void TransponerM(double A[NMAX][NMAX],double AT[NMAX][NMAX],int n, int m)
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++) AT[i][j]=A[j][i];
    }
}
void ProduktMogMT(double A[NMAX][NMAX],double AT[NMAX][NMAX], double P[NMAX][NMAX], int n, int m)
{

    for(int i=0;i<m;i++)
    {
        for (int j=0;j<m;j++)
        {
            P[i][j]=0;
            for(int k=0;k<n;k++) P[i][j] += AT[i][k]*A[k][j];
        }
    }
}
void Totalmatrix(double TM[NMAX][NMAX+1],double A[NMAX][NMAX],double y[NMAX],int m)
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<m;j++) TM[i][j]=A[i][j];
        TM[i][m]=y[i];
    }
}
void Tagsoejle(double A[NMAX][NMAX], double v[NMAX], int n, int k)
{
    for (int i=0;i<n;i++) v[i]=A[i][k];
}

//Vektor:
void RetVektor(double v[NMAX],int n)
{
    int i;
    char Svar;

    cout<<"\nSkal der rettes i vektoren? J/N: ";
    cin>>Svar;

    while(Svar=='J') {
        cout<<"\nIndtast række nr: ";
        cin>>i;i--;
        cout<<"\nIndtast ny værdi: ";
        cin>>v[i];
        cout<<"\nNy vektor givet ved: "<<endl;
        Udskrivv(v,n);
        cout<<"\nSkal der rettes i vektoren igen? J/N: ";
        cin>>Svar;
    }
}
void Indprodukt(double a[NMAX], double b[NMAX], double p[NMAX], int n)
{
    for (int i=0;i<n;i++) p[i] = a[i]*b[i];
}
double Pprodukt(double a[NMAX], double b[NMAX], int n)
{
    double sum = 0;
    for (int i=0;i<n;i++) sum += a[i]*b[i];

    return sum;
}
double Ldiffvek(double a[NMAX], double b[NMAX], double n)
{
    double sum = 0;
    for (int i=0;i<n;i++)
    {
        sum += pow((a[i]-b[i]),2);
    }

    return sqrt(sum);
}
double Lvek(double v[NMAX], int n)
{
    double sum=0;

    for (int i=0;i<n;i++)
    {
        sum += v[i]*v[i];
    }

    return sqrt(sum);
}

//----------------- FUNKTIONER TIL DE ENKELTE METODER -----------------------------------
//Vigtige funktioner:
void BackwardsSubstitution(double TM[NMAX][NMAX+1],double x[NMAX],int n)
{
    double sum;

    x[n-1]=TM[n-1][n]/TM[n-1][n-1];
    for(int i=n-2;i>=0;i--)
    {
        sum=0;
        for(int j=i+1;j<n;j++) sum+=TM[i][j]*x[j];
        x[i]=(TM[i][n]-sum)/TM[i][i];
    }
}

//Metode 1:
void DelvisPivotering(double TM[NMAX][NMAX+1],int j,int n)
{
    double a;
    int r;

    r=j;
    for(int i=j+1;i<n;i++)
    {
        if(fabs(TM[r][j])<fabs(TM[i][j])) r=i;
    }
    if(r!=j){
        for(int k=j;k<=n;k++)
        {
            a=TM[j][k];
            TM[j][k]=TM[r][k];
            TM[r][k]=a;
        }
    }
}
void Gauss(double TM[NMAX][NMAX+1],int m,int &bs)
{
	int i,j,k;
	double factor,eps=0.00000001;
	bs=1;
	for(j=0;j<=m-2;j++)
	{
		DelvisPivotering(TM,j,m);
		if(fabs(TM[j][j])<eps)
		{
			bs=0;
			cout<<"Der forekommer singularitet.";
			break;
		}
		for(i=j+1;i<=m-1;i++)
		{
			factor=-TM[i][j]/TM[j][j];
			TM[i][j]=0;
			for(k=j+1;k<=m;k++)
			{
				TM[i][k]=TM[i][k]+factor*TM[j][k];
			}
		}
	}
	if(fabs(TM[m-1][m-1])<eps)
	{
		bs=0;
		cout<<"Der forekommer singularitet.";
	}
}
void Metode1(double A[NMAX][NMAX],double b[NMAX],double x[NMAX],int n, int m)
{
    double TM[NMAX][NMAX+1], AT[NMAX][NMAX], P[NMAX][NMAX], p[NMAX];
    int bs;

    TransponerM(A,AT,n,m);
    cout<<"\nDen transponeret matrix AT er givet ved: "<<endl;
    UdskrivM(AT,m,n);
    ProduktMogMT(A,AT,P,n,m);
    cout<<"\nMatrix produktet ATA er givet ved: "<<endl;
    UdskrivM(P,m,m);
    ProduktMogv(AT,b,p,m,n);
    cout<<"\nMatrix/vektor produktet ATb er givet ved: "<<endl;
    Udskrivv(p,m);
    Totalmatrix(TM,P,p,m);
    Gauss(TM,m,bs);
    if(bs==1) BackwardsSubstitution(TM,x,m);
    else cout << "Det lineaere ligningssystem kunne ikke bestemmes." << endl;
    cout<<"\nx* givet ved: "<<endl;
    Udskrivv(x,m);
}

//Metode 2:
void Metode2(double A[NMAX][NMAX], double b[NMAX], double Q[NMAX][NMAX], double R[NMAX][NMAX], int n, int m)
{
    double TM[NMAX][NMAX+1], QT[NMAX][NMAX];
    double y[NMAX], x[NMAX];
    TransponerM(Q,QT,n,m);
    ProduktMogv(QT,b,y,m,n);
    Totalmatrix(TM,R,y,m);
    BackwardsSubstitution(TM,x,m);
    cout<<"\nx* givet ved: "<<endl;
    Udskrivv(x,m);
}
void DanQ(double A[NMAX][NMAX], double Q[NMAX][NMAX], double R[NMAX][NMAX], int n, int k)
{
    for (int i=0;i<n;i++)
    {
        Q[i][k] = (1/R[k][k])*A[i][k];
    }
}
void QRmodGS(double A[NMAX][NMAX], double Q[NMAX][NMAX], double R[NMAX][NMAX], int n, int m)
{
    double AA[NMAX][NMAX], v[NMAX];

    for (int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++) AA[i][j] = A[i][j];
    }
    for(int k=0;k<(m-1);k++)
    {
        Tagsoejle(AA,v,n,k);
        R[k][k] = Lvek(v,n);
        DanQ(AA,Q,R,n,k);
        for (int j=k+1;j<m;j++)
        {
            R[k][j] = 0;
            for (int i=0;i<n;i++) R[k][j] += Q[i][k]*AA[i][j];
            for (int i=0;i<n;i++) AA[i][j] = AA[i][j]-(R[k][j]*Q[i][k]);
        }
        Tagsoejle(AA,v,n,(m-1));
        R[m-1][m-1] = Lvek(v,n);
        DanQ(AA,Q,R,n,(m-1));
    }
}

//Metode 3:
void Metode3(double A[NMAX][NMAX], double H[NMAX][NMAX], double b[NMAX], int n, int m)
{
    double grad[NMAX], sk[NMAX], x[NMAX], xgl[NMAX];
    double eps, alph=0;
    int k=-1, N;

    ParameterMetode3(x,eps,N,m);
    do
    {
        k++;
        if (k == 0) DanHesse(H,m);
        else Hesse(A,b,H,x,xgl,n,m);
        Gradient(A,b,x,grad,n,m);
        BFGS(H,grad,sk,m);
        GoldenSection(A,b,x,sk,n,m,alph);
        for (int i=0;i<m;i++){
            xgl[i] = x[i];
            x[i] = x[i] + alph*sk[i];
        }
    } while (Ldiffvek(x,xgl,m)>eps  &&  Lvek(grad,m)>eps &&  k<N);
    cout<<"\nAntal iterationer: "<<k<<endl;
    cout<<"\nVektoren er givet ved:"<<endl;
    Udskrivv(x,m);
}
void ParameterMetode3(double x[NMAX], double &eps,  int &N, int m)
{
    cout << "\nAngiv dit startpunkt, x:" << endl;
    Angivv(x,m);
    cout << "\nAngiv maksimalt antal af iterationer, N = ";
    cin >> N;
    cout << "\nAngiv tolerancen, eps = ";
    cin >> eps;
}
double f(double A[NMAX][NMAX], double b[NMAX], double x[NMAX], int n, int m)
{
    double AT[NMAX][NMAX], P[NMAX][NMAX], p[NMAX];
    double led1 = 0, led2 = 0, led3 = 0;

    for (int i=0;i<n;i++) led1 += b[i]*b[i];
    ProduktMogv(A,x,p,n,m);
    for (int i=0;i<n;i++) led2 += 2*b[i]*p[i];
    TransponerM(A,AT,n,m);
    ProduktMogMT(A,AT,P,n,m);
    ProduktMogv(P,x,p,n,m);
    for (int i=0;i<n;i++) led3 += x[i]*p[i];

    return led1-led2+led3;
}
double falpha(double A[NMAX][NMAX], double b[NMAX], double x[NMAX], double sk[NMAX], double alpha, int n, int m)
{
    double v[NMAX];

    for(int i=0;i<n;i++) v[i]=x[i]+alpha*sk[i];

    return f(A,b,v,n,m);
}
void Gradient(double A[NMAX][NMAX], double b[NMAX], double x[NMAX], double grad[NMAX], int n, int m)
{
    double AT[NMAX][NMAX], P[NMAX][NMAX];
    double p1[NMAX], p2[NMAX];

    TransponerM(A,AT,n,m);
    ProduktMogMT(A,AT,P,n,m);
    ProduktMogv(P,x,p1,n,m);
    ProduktMogv(AT,b,p2,m,n);
    for (int i=0;i<m;i++) grad[i] = 2*p1[i] - 2*p2[i];
}
void DanHesse(double H[NMAX][NMAX], int m)
{
    for (int i=0;i<m;i++)
    {
        for (int j=0;j<m;j++)
        {
            if (i==j) H[i][j] = 1;
            else H[i][j] = 0;
        }
    }
}
void Hesse(double A[NMAX][NMAX], double b[NMAX], double H[NMAX][NMAX], double xny[NMAX], double xgl[NMAX], int n, int m)
{
    double q[NMAX], t[NMAX], g[NMAX], gradgl[NMAX], gradny[NMAX], p[NMAX];

    Gradient(A,b,xgl,gradgl,n,m);
    Gradient(A,b,xny,gradny,n,m);
    for (int i=0;i<m;i++)
    {
        t[i] = xny[i] - xgl[i];
        g[i] = gradny[i] - gradgl[i];
    }
    ProduktMogv(H,g,p,m,m);
    for (int i=0;i<m;i++) q[i] = t[i]-p[i];
    for (int i=0;i<m;i++)
    {
        for (int j=0;j<m;j++) H[i][j] = H[i][j] + (q[i]*t[j]+t[i]*q[j])/Pprodukt(g,t,m) - (Pprodukt(q,g,m)/pow(Pprodukt(g,t,m),2))*t[i]*t[j];
    }
}
void BFGS(double H[NMAX][NMAX], double grad[NMAX], double sk[NMAX], int m)
{
    double v[NMAX];

    for (int i=0;i<m;i++)
    {
        v[i] = 0;
        for (int j=0;j<m;j++) v[i] += -H[i][j]*grad[j];
    }
    for (int i=0;i<m;i++) sk[i] = v[i]/Lvek(v,m);
}
void IndkredsningsAlg(double A[NMAX][NMAX], double vb[NMAX], double x[NMAX], double sk[NMAX], double d, int n, int m, double &a, double &b)
{
    double x1, x2, x3;

    x1=0;
    x2=x1+d;
    while (falpha(A,vb,x,sk,x2,n,m)>=falpha(A,vb,x,sk,x1,n,m))
    {
        d=0.1*d;
        x2=x1+d;
    }
    d=2*d;
    x3=x2+d;
    while (falpha(A,vb,x,sk,x2,n,m)>falpha(A,vb,x,sk,x3,n,m))
    {
        x1=x2;
        x2=x3;
        d=2*d;
        x3=x2+d;
    }
    a=x1;
    b=x3;
}
void GoldenSection(double A[NMAX][NMAX], double b[NMAX], double x[NMAX],double sk[NMAX],int n,int m,double &alph)
{
    double a1, a2, a3, a4, c, eps=0.00000001,d=0.1;

    IndkredsningsAlg(A,b,x,sk,d,n,m,a1,a4);
    c=(sqrt(5.0)-1)/2;
    a2=a1+(1-c)*(a4-a1);
    a3=a4-(1-c)*(a4-a1);
    do
    {
        if (falpha(A,b,x,sk,a2,n,m)>=falpha(A,b,x,sk,a3,n,m))
        {
            a1=a2;
            a2=a3;
            a3=a4-(1-c)*(a4-a1);
        }
        else
        {
            a4=a3;
            a3=a2;
            a2=a1+(1-c)*(a4-a1);
        }
    } while((a4-a1)>=eps);
    alph=0.5*(a4+a1);
}

//Analyse af Metode 3:
void AnalyseMetode3(double R[NMAX][NMAX], double P[NMAX][NMAX], double H[NMAX][NMAX], int m)
{
    double Rinv[NMAX][NMAX];
    DanRinv(R,Rinv,m);
    RinvogRinvT(Rinv,P,m);
    cout << endl;
    cout<<"Matricen dannet fra 1/2 * RInv * RInvt er givet ved: "<<endl;
    UdskrivM(P,m,m);
    cout << endl;
    cout<<"Iterationsmatricen H(m+1) er givet ved: "<<endl;
    UdskrivM(H,m,m);
}
void Danenhedsvek(double e[NMAX], int m, int k)
{
    for (int i=0;i<m;i++)
    {
        if (i == k) e[i] = 1;
        else e[i] = 0;
    }
}
void DanRinv(double R[NMAX][NMAX], double Rinv[NMAX][NMAX], int m)
{
    double TM[NMAX][NMAX+1], x[NMAX], e[NMAX];

    for (int i=0;i<m;i++)
    {
        Danenhedsvek(e,m,i);
        Totalmatrix(TM,R,e,m);
        BackwardsSubstitution(TM,x,m);
        for (int j=0;j<m;j++) Rinv[j][i] = x[j];
    }
}
void RinvogRinvT(double Rinv[NMAX][NMAX], double P[NMAX][NMAX], int m)
{
    double RinvT[NMAX][NMAX];

    TransponerM(Rinv,RinvT,m,m);
    ProduktMogMT(RinvT,Rinv,P,m,m);
    for (int i=0;i<m;i++)
    {
        for (int j=0;j<m;j++) P[i][j] = P[i][j]/2;
    }
}
