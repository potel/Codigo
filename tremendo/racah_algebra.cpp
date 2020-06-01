#include <armadillo>
using namespace arma;
using namespace std;
#include "tremendo.h"
#include "structs.h"
#include "definiciones.h"

double fact(int n)
{
  if (n>0) return log(n)+fact(n-1);
  else return 0.;
}

inline double phase(int i)
{
  return pow(-1.,i);
}
inline int frac(double x)
{
  return abs(x-floor(x))>1.e-15;
}
inline int fail3(double x,double y,double z)
{
  return ((frac(x+y+z)) || (x>(y+z)) || (x<abs(y-z)));
}
double wig6j(double a,double b,double c, double d, double e, double f)
{
  double ic,id,ie,ig,ia,ib,ih,m,mm,
    mup,t,n,s,it,iu,iv,iw,ta,tb,xd;
  if (fail3(a,b,c)) return 0.;
  if (fail3(a,e,f)) return 0.;
  if (fail3(b,d,f)) return 0.;
  if (fail3(c,d,e)) return 0.;
  ic=a+b-c;
  id=e+d-c;
  ie=a+e-f;
  ig=b+d-f;
  ib=c+f-a-d;
  ia=c+f-b-e;
  ih=a+b+e+d+1;
  mm=min(min(ih,ic),min(id,ie));
  m=min(mm,ig);
  if (m<0.) return 0.;
  mup=min(min(ia,ib),0.);
  t=phase(m);
  n=m;
  s=t;
  m=m-1;
  while (m+mup>=0.)
    {
      ta=(ia+m+1.)*(ib+m+1.)*(ih-m)*(m+1.);
      tb=(ic-m)*(id-m)*(ie-m)*(ig-m);
      t=-t*ta/tb;
      s=s+t;
      m=m-1;
    }
  it=a+b+c+1.;
  iu=a+e+f+1.;
  iv=b+d+f+1.;
  iw=c+d+e+1.;
  xd=0.5*(fact(ic)+fact(ie+ib)+fact(ia+ig)+fact(ie)+fact(ib+ic)+fact(ia+id)+fact(ig)+fact(ic+ia)
          +fact(id+ib)+fact(id)+fact(ia+ie)+fact(ib+ig)-fact(it)-fact(iu)-fact(iv)-fact(iw))
    +fact(ih-n)-fact(n)-fact(ia+n)-fact(ib+n)-fact(ic-n)-fact(id-n)-fact(ie-n)-fact(ig-n);
  return phase(ih-1.)*s*exp(xd);
}

double wig9j(double a,double b,double c, double d, double e, double f, double g, double h, double z)
{
  double s,x,xa,xb,xc,k;
  if (fail3(a,b,c)) return 0.;
  if (fail3(d,e,f)) return 0.;
  if (fail3(c,f,z)) return 0.;
  if (fail3(a,d,g)) return 0.;
  if (fail3(b,e,h)) return 0.;
  if (fail3(g,h,z)) return 0.;
  s=0.;
  xa=abs(a-z);
  xb=abs(d-h);
  xc=abs(b-f);
  x=xa;
  if ((x-xb)<0.)  x=xb;
  else
    {
      //cout<<"quillo 1\n";
      if ((x-xc)<0.) x=xc;
      //else
        {
          //  cout<<"quillo 2\n";
          if((x-a-z)<=0.)
            {
              //  cout<<"quillo 3\n";
              if((x-d-h)<=0.)
                {
                  //  cout<<"quillo 4\n";
                  if((x-b-f)<=0) s+=(2.*x+1.)*wig6j(a,z,x,h,d,g)*wig6j(b,f,x,d,h,e)*wig6j(a,z,x,f,b,c);
                  //cout<<x<<"  "<<s<<"\n";
                  x++;
                }
            }
          else
            if (s==0.) return 0.;
            else k=2.*(a+b+d+f+h+z);
        }
    }
  return phase(k)*s;
}


















