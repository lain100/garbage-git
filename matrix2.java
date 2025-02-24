public class foo {
  public static void main(String[] args) {
    double[][] x0 = {
          {0},{0}};
  System.out.println(
      simpson(
        -1, 2, 1000));
    double[][] a = {
            {1, -2, -4},
            {2, -4, 3},
            {3, 4, -5}};
    double[][] b = {
            {-4},
            {3},
            {5}};
            
    // showMatrix(
    //     gauss( a, b ));
  }
  
  public static void bisection(double a, double b) {
    if ( f(a) * f(b) > 0 ) {
      System.out.println("No solution exists for the given range of initial values.");
      System.exit(1);
    }
    
    double m = 0.0, e = 1e-8;
    while (true) {
      m = (a + b) / 2.0;
      if ( Math.abs( f(m) ) < e || Math.abs(b - a) < e) break;
      if (f(a) * f(m) > 0)  a = m;
      else          b = m;
    }
    System.out.printf("x = %f f(x) = %f\n", m, f(m) );
  }
  
  public static void falsePosition(double a, double b) {
    if ( f(a) * f(b) > 0 ) {
      System.out.println("No solution exists for the given range of initial values.");
      System.exit(1);
    }
    
    double m = 0.0, e = 1e-8;
    while (true) {
      m = (a * f(b) - b * f(a)) / (f(b) - f(a));
      if ( Math.abs( f(m) ) < e || Math.abs(b - a) < e) break;
      if (f(a) * f(m) > 0)  a = m;
      else          b = m;
    }
    System.out.printf("x = %f f(x) = %f\n", m, f(m));
  }
  
  public static void newton(double x0) {
    double e = 1e-8;
    double h = 1e-8;
    int count = 0;
    int max = 100;
    
    double x = x0;
    for (; count < max; count++) {
      if ( df(x, h) == 0.0 ) {
        System.out.println("Cannot be calculated because df(x) = 0.");
        System.exit(1);
      }
      x -= f(x) / df(x, h);
      if ( Math.abs(f(x)) < e ) break;
    }
    
    if (count == max) {
      System.out.println("The solution did not converge.");
    }
    
    System.out.printf("x = %f f(x) = %f\n", x, f(x));
  }
  
  public static void newton(double[][] x0) {
    double e = 1e-8;
    double h = 1e-8;
    int count = 0;
    int max = 100;
    
    double[][] x = copyMatrix(x0);
    for (; count < max; count++) {
      double[][] f = f(x);
      double[][] df = df(x, h);
      double[][] y = gauss(df, f);
      double[][] xi = subMatrix(x, y);
      
      if (maxNorm(x, xi) < e) break;
      x = copyMatrix(xi);
    }
    
    if (count == max) {
      System.out.println("The solution did not converge.");
    }
    
    showMatrix(x);
    showMatrix(f(x));
  }
  
  public static double maxNorm(double[][] x, double[][] xi) {
    int m = x.length;
    int n = x[0].length;
    double norm = 0.0;
    
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        double tmp = Math.abs(x[i][j] - xi[i][j]);
        if (norm < tmp) norm = tmp;
      }
    }
    return norm;
  }
  
  public static void showMatrix(double[][] a) {
    int m = a.length;
    int n = a[0].length;
    
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        System.out.print(a[i][j] + " ");
      }
      System.out.println();
    }
    System.out.println();
  }
  
  public static double[][] subMatrix(double[][] a, double[][] b) {
    int m = a.length;
    int n = a[0].length;
    double[][] res = new double[m][n];
    
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        res[i][j] = a[i][j] - b[i][j];
      }
    }
    
    return res;
  }
  
  public static double[][] NaN_Matrix(double[][] x) {
    int m = x.length;
    int n = x[0].length;
    
    for (int i = 0; i < m; i++) {      
      for (int j = 0; j < n; j++) {
        x[i][j] = Double.NaN;
      }
    }
    
    return x;
  }
  
  public static double[][] copyMatrix(double[][] a) {
    int m = a.length;
    int n = a[0].length;
    double[][] b = new double[m][n];
    
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        b[i][j] = a[i][j];
      }
    }
    return b;
  }
  
  public static double[][] gauss(double[][] _a, double[][] _b) {
    int n = _a.length;
    double[][] a = copyMatrix( _a );
    double[][] x = copyMatrix( _b );
    
    for (int i = 0; i < n; i++) {  
      double max = Math.abs(a[i][i]);
      int index = i;
      
      for (int j = i + 1; j < n; j++) {
        double val = Math.abs(a[j][i]);
        if (max < val) {
          index = j;
          max = val;
        }
      }
      
      if (max < 1e-24) {
        return NaN_Matrix( x );
      }
      
      if (index != i) {
        for (int j = i; j < n; j++) {
          double temp = a[index][j];
          a[index][j] = a[i][j];
          a[i][j] = temp;
        }
        
        double temp = x[index][0];
        x[index][0] = x[i][0];
        x[i][0] = temp;
      }
      
      for (int j = i + 1; j < n; j++) {
        double ratio = a[j][i] / a[i][i];
        for (int k = i; k < n; k++) {
          a[j][k] -= ratio * a[i][k];
        }
        x[j][0] -= ratio * x[i][0];
      }
    }
    
    for (int i = n - 1; i >= 0; i--) {      
      for (int j = i + 1; j < n; j++) {
        x[i][0] -= a[i][j] * x[j][0];
      }
      x[i][0] /= a[i][i];
    }
    
    return x;
  }
  
  public static double[][] f(double[][] x) {
    double[][] f = new double[2][1];
    
    f[0][0] = x[0][0]*(x[0][0] - 3)*(x[0][0] - 3) - x[1][0]*x[1][0] + 1;
    f[1][0] = x[1][0]*(Math.exp(x[0][0]) - Math.exp(-x[0][0])) - 4*x[1][0]*x[1][0] + 1;
    
    return f;
  }
  
  public static double[][] df(double[][] x, double h) {
    int n = x.length;
    double[][] df = new double[n][n];
    
    for (int j = 0; j < n; j++) {
      double[][] xh = copyMatrix(x);
      xh[j][0] += h;
      double[][] f = f(x);
      double[][] fh = f(xh);
        
      for (int i = 0; i < n; i++) {
        df[i][j] = (fh[i][0] - f[i][0]) / h;
      }
    }
    
    return df;
  }
  
  public static double f(double x) {
    return -3*x*x+4*x+2;
  }
  
  public static double df(double x, double h) {
    return (f(x + h) - f(x)) / h;
  }
  
  public static double trapezoid(double a, double b, int n) {
    double h = (b - a) / n;
    double x = a;
    double tmp = 0.0;
    
    for (int i = 1; i <= n - 1; i++) {
      x += h;
      tmp += f(x);
    }
    
    double integral = h * ((f(a) + f(b)) / 2.0 + tmp);
    
    return integral;
  }
  
  public static double midpoint(double a, double b, int n) {
    double h = (b - a) / n;
    double x = a + h / 2.0;
    double integral = 0.0;
    
    for (int i = 0; i <= n - 1; i++) {
      integral += h * f(x);
      x += h;
    }

    return integral;
  }
  
  public static double simpson(double a, double b, int n) {
    double h = (b - a) / n;
    double x = a;
    double tmp1 = 0.0, tmp2 = 0.0;
    
    for (int i = 1; i <= n/2 - 1; i++) {
      x += 2*h;
      tmp1 += f(x - h);
      tmp2 += f(x);
    }
    
    tmp1 += f(b - h);
    
    double integral = h/3.0 * (f(a) + 4.0*tmp1 + 2.0*tmp2 + f(b));
    
    return integral;
  }
}
