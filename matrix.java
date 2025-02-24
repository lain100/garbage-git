public class matrix {
  public static void main(String[] args) {
    double[][] a = {
      {1, 2, -3},
      {-2, -3, 4},
      {3, -4, -1}};
    double[][] b = {
      {-6},
      {5},
      {4}};
     showMatrix(lu( a, b ));
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
  
  public static double det(double[][] _a) {
    int n = _a.length;
    double detA = 1.0;
    double[][] a = copyMatrix( _a );
    
    for (int i = 0; i < n; i++) {  
      double max = Math.abs(a[i][i]);
      int index = i;
      
      for (int j = i + 1; j < n; j++) {
        double x = Math.abs(a[j][i]);
        if (max < x) {
          index = j;
          max = x;
        }
      }
      
      if (max < 1e-24) return 0.0;
      
      if (index != i) {
        detA *= -1;
        for (int j = i; j < n; j++) {
          double temp = a[index][j];
          a[index][j] = a[i][j];
          a[i][j] = temp;
        }
      }
      
      for (int j = i + 1; j < n; j++) {
        double ratio = a[j][i] / a[i][i];
        for (int k = i; k < n; k++) {
          a[j][k] -= ratio * a[i][k];
        }
      }
    }
    
    for (int i = 0; i < n; i++) {
      detA *= a[i][i];
    }
    return detA;
  }
  
  public static double[][] cramer(double[][] a, double[][] b) {
    int n = a.length;
    double detA = det( a );
    double[][] x = new double[n][1];
    
    if (detA == 0.0) return NaN_Matrix( x );
    
    for (int i = 0; i < n; i++) {
      double[][] ai = copyMatrix( a );
      for (int j = 0; j < n; j++) {
        ai[j][i] = b[j][0];
      }
      x[i][0] = det( ai ) / detA;
    }
    
    return x;
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
  
  public static double[][] lu(double[][] a, double[][] b) {
    int n = a.length;
    double[][] lu = luDecomp( a );
    double[][] x = copyMatrix( b );
    
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        x[i][0] -= lu[i][j] * x[j][0];
      }
    }
    
    for (int i = n - 1; i >= 0; i--) {
      for (int j = i + 1; j < n; j++) {
        x[i][0] -= lu[i][j] * x[j][0];
      }
      x[i][0] /= lu[i][i];
    }
    
    return x;
  }
  
  public static double[][] luDecomp(double[][] _a) {
    int n = _a.length;
    double[][] a = copyMatrix( _a );
    
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
          a[j][i] /= a[i][i];
        for (int k = i + 1; k < n; k++) {
          a[j][k] -= a[i][k] * a[j][i];
        }
      }
    }
    
    return a;
  }
  
  public static double[][] inv(double[][] _a) {
    int n = _a.length;
    double[][] a = copyMatrix( _a );
    double[][] lu = luDecomp( a );
    double[][] inv = new double[n][n];
    
    for (int k = 0; k < n; k++) {
      for (int i = 0; i < n; i++) {
        inv[i][k] = 0.0;
      }
      inv[k][k] = 1.0;
      
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
          inv[i][k] -= lu[i][j] * inv[j][k];
        }
      }
      
      for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
          inv[i][k] -= lu[i][j] * inv[j][k];
        }
        inv[i][k] /= lu[i][i];
      }
    }
    
    return inv;
  }
}
