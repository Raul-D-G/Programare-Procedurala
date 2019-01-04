#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

//Tip de date utilizat pentru elimiarea non-maximelor
typedef struct
{
    double corr;
    int x, y, latime, inaltime;
    int R, G, B;
} Detectie;

//Generatorul de numere aleatoare
uint32_t xorshift32(uint32_t *seed)
{
    uint32_t x = (*seed);
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    (*seed) = x;
    return x;
}
//Interschimbare elemente
void changeValues (int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}
//Obtinere Permutare aleatoare
void shuffleRandon ( int **arr1, int n , uint32_t *R)
{
    int i;
    for ( i = n-1; i > 0; i--)
    {

        int j = R[n-i] % (i+1);
        changeValues(&(*arr1)[i], &(*arr1)[j]);
    }
}
//Incarca o imagine in memoria interna
unsigned char * incarcare_imagine(char * nume_imagine, int **L, int *W, int *H)
{
    FILE *f = fopen(nume_imagine,"rb");
    if(f == NULL)
    {
        printf("nu am gasit imaginea sursa din care citesc");
        return 0;
    }

    int h, i, j;
    unsigned char *header;
    unsigned int latime_img, inaltime_img;

    fseek(f, 0, SEEK_SET);

    header=(unsigned char *)calloc(54,sizeof(char));
    for(h = 0; h < 54; h++)
    {
        fread(&header[h], 1, 1, f);
    }

    fseek(f, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, f);
    fread(&inaltime_img, sizeof(unsigned int), 1, f);

    int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;


    (*L)=(int *)calloc(latime_img*inaltime_img, sizeof(int));

    (*H)=inaltime_img;
    (*W)=latime_img;

    fseek(f, 54, SEEK_SET);

    unsigned char pixel[3];
    int aux;
    for(i = inaltime_img-1; i >= 0; i--)
    {
        for(j = 0; j < latime_img; j++)
        {
            fread(pixel, 3, 1, f);
            aux = pixel[2];
            aux = aux<<8;
            aux ^= pixel[1];
            aux = aux<<8;
            aux ^= pixel[0];
            (*L)[(i * latime_img) + j] = aux;
        }
        fseek(f, padding, SEEK_CUR);
    }


    return header;
    fclose(f);
}
//Salveaza o imagine din memoria interna in memoria externa
void salvare_imagine(char *destinatie,int *L,int latime_img, int inaltime_img,  unsigned char *header)
{
    int i,j;

    FILE *g=fopen(destinatie, "wb");

    int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    for(i=0; i<54; i++)
    {
        fwrite(&header[i],1,1,g);
        fflush(g);
    }

    int h=0;


    for(i = inaltime_img-1; i >= 0; i--)
    {
        for(j = 0; j < latime_img; j++)
        {
            fwrite(&L[(i*latime_img)+j],3,1,g);
            fflush(g);
            h++;

            if(h == latime_img)
            {
                char O=0;
                int t;
                for(t=0; t<padding; t++)
                {
                    fwrite(&O, 1, 1, g);
                    fflush(g);
                }
                h=0;
            }
            fflush(g);
        }
    }

    free(L);
    fclose(g);
}
//Operatia XOR intre 2 pixeli, primul se modifica
void XOR_pixel(int *x,int y)
{
    unsigned char p1[3], p2[3], i, octet;
    int aux;
    for(i = 0; i < 3; i++)
    {
        octet = (*x);
        (*x) = (*x)>>8;
        p1[i] = octet;
        octet = y;
        y = y>>8;
        p2[i] = octet;
    }

    for(i = 0; i < 3; i++)
        p1[i] = p1[i] ^ p2[i];

    aux = p1[2];
    aux = aux<<8;
    aux ^= p1[1];
    aux = aux<<8;
    aux ^= p1[0];
    (*x) = aux;
}
//Cripteaza o imagine "sursa" si o salveaza in "destinatie" pe baza informatiilor din "cheie_secreta"
void criptare(char *sursa, char *destinatie, char *cheia_secreta)
{
    int i, inaltime_img, latime_img;
    unsigned int SV;
    uint32_t *R=NULL, R0;
    int *L=NULL;
    int *Per;
    unsigned char *header=NULL;

    FILE *f=fopen(sursa, "rb");
    if(f == NULL)
    {
        printf("nu am gasit imaginea sursa din care citesc");
        return ;
    }

    FILE *g=fopen(destinatie, "wb");
    if(g == NULL)
    {
        printf("eroare la creearea imaginii");
        return ;
    }

    FILE *K=fopen(cheia_secreta, "r");
    if(K == NULL)
    {
        printf("nu am gasit cheia secreta");
        return ;
    }

    fscanf(K, "%u %u", &R0, &SV);

    header=incarcare_imagine(sursa, &L, &latime_img, &inaltime_img);

    R=(unsigned int *)calloc((2*latime_img*inaltime_img)-1,sizeof(unsigned int));


    //XORSHIFT32
    for(i=1; i<=2*latime_img*inaltime_img-1; i++)
        R[i]=xorshift32(&R0);

    //Obtinere permutare
    Per=(int *)calloc(latime_img*inaltime_img,sizeof(int));
    for(i=0; i<latime_img*inaltime_img; i++)
        Per[i]=i;

    shuffleRandon(&Per, latime_img*inaltime_img, R);


    //Imagine cu pixelii permutati

    int *Ll;
    Ll=(int *)calloc(latime_img*inaltime_img,sizeof(int));
    for(i=0; i<latime_img*inaltime_img; i++)
        Ll[Per[i]]=L[i];

    for(i=0; i<latime_img*inaltime_img; i++)
        L[i]=Ll[i];
    free(Ll);



    //Criptare
    for(i=0; i<latime_img*inaltime_img; i++)
    {
        if(i==0)
        {
            XOR_pixel(&L[i],SV);
            XOR_pixel(&L[i],R[latime_img*inaltime_img+i]);
        }
        else
        {
            XOR_pixel(&L[i],L[i-1]);
            XOR_pixel(&L[i],R[latime_img*inaltime_img+i]);
        }
    }

    salvare_imagine(destinatie, L, latime_img, inaltime_img,header);
    printf("Criptare realizata!\n");

    fclose(f);
    fclose(g);
    fclose(K);
}
//Decripteaza o imagine "imagine_criptata" si o salveaza in "imagine_initiala" pe baza informatiilor din "cheie_secreta"
void Decriptare(char *imagine_initiala, char *imagine_criptata, char *cheia_secreta)
{
    int i, inaltime_img, latime_img;
    unsigned int SV;
    uint32_t *R=NULL, R0;
    int *L=NULL;
    int *Per,*Per_invers;
    unsigned char *header=NULL;

    FILE *K=fopen(cheia_secreta, "r");
    if(K == NULL)
    {
        printf("nu am gasit cheia secreta");
        return ;
    }

    fscanf(K, "%u %u", &R0, &SV);


    header=incarcare_imagine(imagine_criptata, &L, &latime_img, &inaltime_img);

    R=(unsigned int *)calloc((2*latime_img*inaltime_img)-1,sizeof(unsigned int));


    //XORSHIFT32
    for(i=1; i <= 2*latime_img*inaltime_img-1; i++)
        R[i]=xorshift32(&R0);


    //Obtinere permutare
    Per=(int *)calloc(latime_img*inaltime_img,sizeof(int));
    for(i=0; i<latime_img*inaltime_img; i++)
        Per[i]=i;

    shuffleRandon(&Per, latime_img*inaltime_img, R);

    Per_invers=(int *)calloc(latime_img*inaltime_img,sizeof(int));

    for(i=0; i < latime_img*inaltime_img; i++)
        Per_invers[Per[i]]=i;
    free(Per);


    int *Lk=(int *)calloc(latime_img*inaltime_img,sizeof(int));
    for(i = 0; i < latime_img*inaltime_img; i++)
        Lk[i]=L[i];

    XOR_pixel(&L[0], SV);
    XOR_pixel(&L[0], R[latime_img*inaltime_img]);
    for(i = 1; i < latime_img*inaltime_img; i++)
    {
        XOR_pixel(&L[i], Lk[i-1]);
        XOR_pixel(&L[i], R[latime_img*inaltime_img+i]);
    }

    free(Lk);


//            Imagine cu pixelii permutati
    int *Ll;
    Ll=(int *)calloc(latime_img*inaltime_img,sizeof(int));
    for(i = 0; i < latime_img*inaltime_img; i++)
        Ll[Per_invers[i]] = L[i];

    for(i = 0; i < latime_img*inaltime_img; i++)
        L[i] = Ll[i];
    free(Ll);


    salvare_imagine(imagine_initiala,L,latime_img,inaltime_img,header);
    printf("Decriptare realizata!\n");

    fclose(K);
}
//Afiseaza pe ecran rezultatele testului Hi^2 pe o imagine:"imagine_testata"
void Hi_2(char *imagine_testata)
{
    unsigned char *header=NULL;
    int *L=NULL;
    int latime_img,inaltime_img,i;
    int *R,*G,*B;
    double F=0.0, Hi_r=0.0, Hi_g=0.0, Hi_b=0.0;

    header=incarcare_imagine(imagine_testata,&L,&latime_img,&inaltime_img);

    free(header);

    R=(int *)calloc(256,sizeof(int));
    G=(int *)calloc(256,sizeof(int));
    B=(int *)calloc(256,sizeof(int));

    unsigned char x;

    for(i = 0; i < latime_img*inaltime_img; i++)
    {
        x = L[i];
        L[i] = L[i]>>8;
        B[x]++;

        x = L[i];
        L[i] = L[i]>>8;
        G[x]++;

        x = L[i];
        R[x]++;
    }

    F = (latime_img*inaltime_img)/256;

    for(i = 0; i < 256; i++)
    {
        Hi_r = Hi_r + (R[i]-F)*(R[i]-F)/F;
        Hi_g = Hi_g + (G[i]-F)*(G[i]-F)/F;
        Hi_b = Hi_b + (B[i]-F)*(B[i]-F)/F;
    }

    printf("Hi-squared test on RGB channels for '%s' :\n",imagine_testata);
    printf("R: %.2f \nG: %.2f \nB: %.2f\n",Hi_r,Hi_g,Hi_b);

    free(L);
    free(header);
    free(R);
    free(G);
    free(B);
}
//Transforma o imagine color :"nume_fisier_sursa" intr-o imagine gri :"nume_fisier_destinatie"
void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie)
{
    FILE *fin, *fout;
    unsigned int dim_img, latime_img, inaltime_img;
    unsigned char pRGB[3], aux;

    //printf("nume_fisier_sursa = %s \n",nume_fisier_sursa);

    fin = fopen(nume_fisier_sursa, "rb");
    if(fin == NULL)
    {
        printf("nu am gasit imaginea sursa din care citesc");
        return;
    }

    fout = fopen(nume_fisier_destinatie, "wb+");

    fseek(fin, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fin);
    //printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);
    //printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n",latime_img, inaltime_img);

    //copiaza octet cu octet imaginea initiala in cea noua
    fseek(fin,0,SEEK_SET);
    unsigned char c;
    while(fread(&c,1,1,fin)==1)
    {
        fwrite(&c,1,1,fout);
        fflush(fout);
    }
    fclose(fin);

    //calculam padding-ul pentru o linie
    int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    //printf("padding = %d \n",padding);

    fseek(fout, 54, SEEK_SET);
    int i,j;
    for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++)
        {
            //citesc culorile pixelului
            fread(pRGB, 3, 1, fout);
            //fac conversia in pixel gri
            aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
            pRGB[0] = pRGB[1] = pRGB[2] = aux;
            fseek(fout, -3, SEEK_CUR);
            fwrite(pRGB, 3, 1, fout);
            fflush(fout);
        }
        fseek(fout,padding,SEEK_CUR);
    }
    fclose(fout);
}
//Incarca o fereastra de anumite dimensiuni dintr-o imagine in memoria interna
void incarcare_fereastra(int x, int y, int *L, int latime_S, int inaltime_S, int latime_img, int **F)
{
    int i, j, t=0;

    (*F)=(int *)calloc(latime_S*inaltime_S,sizeof(int));

    for(i = 0; i < inaltime_S; i++)
        for(j = 0; j < latime_S; j++)
        {
            (*F)[t]  =L[(x+i) * latime_img+(j+y)];
            t++;
        }
}
//Calculeaza corelatia dintre o fereastra si un sablon
double corelatie(int *S, int latime_S, int inaltime_S, int *F)
{
    double rez,n;
    int i;
    unsigned char S_intensitate_pixel, F_intensitate_pixel;
    double S_intensitate_medie, F_intensitate_medie, D_standard_S, D_standard_F;
    double suma_S=0.0, medie_S=0.0,medie_F=0.0, suma_F=0.0, suma=0.0;

    n=latime_S*inaltime_S;

    for(i = 0; i < n; i++)
    {
        S_intensitate_pixel = S[i];
        F_intensitate_pixel = F[i];

        medie_S = medie_S + S_intensitate_pixel;
        medie_F = medie_F + F_intensitate_pixel;
    }
    S_intensitate_medie = medie_S/n;
    F_intensitate_medie = medie_F/n;



    for(i = 0; i < n; i++)
    {
        S_intensitate_pixel = S[i];
        F_intensitate_pixel = F[i];
        suma_S = suma_S + (S_intensitate_pixel - S_intensitate_medie) * (S_intensitate_pixel - S_intensitate_medie);
        suma_F = suma_F + (F_intensitate_pixel - F_intensitate_medie) * (F_intensitate_pixel - F_intensitate_medie);
    }

    double a;
    a=1/(n-1);
    D_standard_S = sqrt(a*suma_S);
    D_standard_F = sqrt(a*suma_F);

    a=(1 / (D_standard_F*D_standard_S));

    for(i = 0; i < n; i++)
    {
        S_intensitate_pixel = S[i];
        F_intensitate_pixel = F[i];
        suma = suma + a * (S_intensitate_pixel - S_intensitate_medie) * (F_intensitate_pixel - F_intensitate_medie);
    }

    rez=(1/n)*suma;
    return rez;
}
//Coloreaza o detectie intr-o imagine incarcata in memoria interna
void colorare_detectie(int **L, int latime_img , Detectie D)
{
    int i, j, x, y,aux;

    x = D.x;
    y = D.y;

    aux = D.R;
    aux = aux<<8;
    aux ^= D.G;
    aux = aux<<8;
    aux ^= D.B;

    for(i = 0; i < D.inaltime; i++)
    {
        (*L)[(x+i)*latime_img+y] = aux;
        (*L)[(x+i)*latime_img+y+D.latime] = aux;
    }

    for(j = 0; j < D.latime; j++)
    {
        (*L)[x*latime_img+y+j] = aux;
        (*L)[(x+D.inaltime)*latime_img+y+j] = aux;
    }
}
//Ruleaza operatia de template matching pe o imagine :"imagine" cu ajutorul unui sablon :"sablon"
Detectie * Template_matching(char *imagine, char *sablon, double Ps,int *t)
{
    int *L = NULL, *S = NULL;
    int *F = NULL;
    double corr;

    Detectie * D;

    int latime_img, inaltime_img, i, j,latime_S,inaltime_S;

    unsigned char *header=NULL, *header_S=NULL;


    grayscale_image(sablon,"sablon_gri.bmp");

    header=incarcare_imagine("imagine_gri.bmp", &L, &latime_img, &inaltime_img);

    header_S=incarcare_imagine("sablon_gri.bmp", &S, &latime_S, &inaltime_S);



    D=(Detectie *)calloc(latime_img*inaltime_img,sizeof(Detectie));


    for(i = 0; i <= inaltime_img-inaltime_S; i++)
        for(j = 0; j <= latime_img-latime_S; j++)
        {
            incarcare_fereastra(i, j, L, latime_S, inaltime_S, latime_img,&F);
            corr=corelatie(S, latime_S, inaltime_S, F);
            if(corr > Ps)
            {
                D[(*t)].corr = corr;
                D[(*t)].x = i;
                D[(*t)].y = j;
                D[(*t)].inaltime=inaltime_S;
                D[(*t)].latime=latime_S;
                (*t)++;
            }
        }



    free(header);
    free(header_S);

    return D;
}
//Functia comparator pentru qsort
int cmp(const void *a,const void *b)
{
    Detectie va = *(Detectie *)a;
    Detectie vb = *(Detectie *)b;
    if(va.corr<vb.corr)
        return 1;
    if(va.corr>vb.corr)
        return -1;
    return 0;
}
//calculeaza aria intersectiei a doua detectii
double intersectie_arie(Detectie a, Detectie b)
{
    int h=a.inaltime,w=a.latime;
    double arie=0.0;

    if(a.x > b.x && a.y >= b.y)
        arie=(b.y + w - a.y)*(b.x + h - a.x);

    if(a.x >= b.x && a.y < b.y)
        arie=(b.x + h - a.x)*(a.y + w - b.y);

    if(a.x < b.x && a.y <= b.y)
        arie=(a.x + h - b.x)*(a.y + w - b.y);

    if(a.x <= b.x && a.y >b.y)
        arie=(a.x + h - b.x)*(b.y + w - a.y);

    return arie;
}
//Calculeaza suprapunerea spatiala a 2 detectii
double suprapunere_arie(Detectie a, Detectie b)
{
    double suprapunere = 0.0;
    int w, h;
    w = a.latime;
    h = a.inaltime;
    suprapunere = intersectie_arie(a,b) / (2*w*h-intersectie_arie(a,b));
    return suprapunere;
}
//Verifica daca 2 detectii se suprapun
int intersectie(Detectie a, Detectie b)
{
    int ok = 1;
    int h = a.inaltime, w = a.latime;
    int x, y;
    x = a.x-b.x;
    y = a.y-b.y;

    if(x < 0)
        x = -x;
    if(y < 0)
        y = -y;
    if(x > h || y > w)
        ok = 0;
    return ok;
}
//Elimina non-maximele din tabloul de detectii si coloreaza detectiile ramase
void eliminare_max_colorare_Detectii(char *imagine, char *nume_sabloane)
{
    int n, i, t, x, j, R, G, B;
    char *nume_sablon;
    FILE *f = fopen(nume_sabloane,"r");
    nume_sablon = (char *)malloc(20*sizeof(char));

    Detectie *D, *M;

    M = (Detectie *)calloc(400*500,sizeof(Detectie));

    grayscale_image(imagine,"imagine_gri.bmp");

    t=0;
    for(x = 0; x < 10; x++)
    {
        fscanf(f, "%s", nume_sablon);
        fscanf(f, "%d %d %d", &R, &G, &B);
        n = 0;
        D = Template_matching(imagine, nume_sablon, 0.5, &n);

        for(i = 0; i < n; i++)
        {
            M[t] = D[i];
            M[t].R = R;
            M[t].G = G;
            M[t].B = B;
            t++;
        }
    }
    qsort(M, t, sizeof(M[0]), cmp);

    for(i = 0; i < t-1; i++)
        for(j = i+1; j < t; j++)
            if(intersectie(M[i], M[j]) == 1)
                if(suprapunere_arie(M[i], M[j]) > 0.2)
                {
                    for(x = j; x < t; x++)
                        M[x] = M[x+1];
                    t--;
                }

    int *L = NULL, W, H;
    unsigned char *header = NULL;
    header = incarcare_imagine("imagine_gri.bmp", &L, &W, &H);

    for(i = 0; i < t; i++)
        colorare_detectie(&L, W, M[i]);

    salvare_imagine("Rezultat.bmp", L, W, H, header);

    printf("Numele imaginii obtinuta in urma operatiei de Template matching este : Rezultat.bmp\n");

    fclose(f);
}

int main()
{
    char img_pentru_criptat[20], img_criptata[20], secret_key_criptare[20], img_pentru_decriptare[20],img_decriptata[20],secret_key_decriptare[20];

    printf("Numele imaginii propusa criptarii: ");
    scanf("%s",img_pentru_criptat);

    printf("Numele imaginii criptate : ");
    scanf("%s",img_criptata);

    printf("Cheia secreta pentru criptare:");
    scanf("%s",secret_key_criptare);

    criptare(img_pentru_criptat,img_criptata,secret_key_criptare);

    printf("Numele imaginii propusa decriptarii : ");
    scanf("%s",img_pentru_decriptare);

    printf("Numele imaginii decriptate : ");
    scanf("%s",img_decriptata);

    printf("Cheia secreta pentru decriptare:");
    scanf("%s",secret_key_decriptare);

    Decriptare(img_decriptata,img_pentru_decriptare,secret_key_decriptare);


    Hi_2(img_decriptata);
    Hi_2(img_criptata);

    char imagine[20];

    printf("Nume imagine pentru template matching :");
    scanf("%s",imagine);

    eliminare_max_colorare_Detectii(imagine,"Nume_sablon.txt");


    return 0;
}
