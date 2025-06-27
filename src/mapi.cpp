#include <iostream>
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <math.h>

#if defined(__GNUC__)
#  pragma GCC optimize ("O3")
#endif

// See: s2-earth.R for using same earth radius
#define EARTHRADIUS 6371010.0L

/*
NOTE: "using namespace" are commented out and replaced by the syntax namespace::object in order to avoid ambiguities between std and Rcpp.
WARNING: DO NOT UNCOMMENT!
using namespace std;
using namespace Rcpp;
*/


//' @noRd
// //' Number of intersections between ellipses and grid cells
// //' 
// //' From the list of integer vectors which represents the intersection of grid cells and ellipses, this function counts the total number of matches (items) in this list of vectors.
// //'
// //' @param inter A list of integer vector containing intersections between cells and ellipses
// //' @return the total number of items (sum of lengths of all vectors in the list).
// [[Rcpp::export(.countMatches_cpp)]]
long countMatches_cpp(Rcpp::List inter) {
  unsigned int n = inter.size();
  long resu = 0L;
  //Rcout << "n=" << n << "\n";
  for (unsigned int i=0; i<n; i++) {
    Rcpp::List ells = inter(i);
    long ne = ells.size();
    resu += ne;
  }
  return(resu);
}

//' @noRd
// //' Retrieve weight and values for each cell from overlapping ellipses
// //'
// //' From the list of intersections, cells and ellipses' features, build the list of (value,weight) pairs for each cell
// //'
// //' @param cells the vector of cell
// //' @param inter the vector of intersections cell,ellipses
// //' @param weights the vector of ellipses weights
// //' @param values the vector of ellipses values
// //' @return a list of results by cell which contains lists of (value,weight) pairs 
// [[Rcpp::export(.getValues_cpp)]]
Rcpp::List getValues_cpp(Rcpp::NumericVector cells, Rcpp::List inter, Rcpp::DoubleVector weights, Rcpp::DoubleVector values) {
  unsigned int n = inter.size();
  Rcpp::List resu(0);
  for (unsigned int i=0; i<n; i++) {
    Rcpp::List ells = inter(i);
    unsigned int ne = ells.size();
    if (ne == 0) {
        Rcpp::NumericMatrix elem(0,0);
        resu.push_back(elem);
    } else {
      Rcpp::NumericMatrix elem(ne,2);
      for (unsigned int j=0; j<ne; j++) {
        int ie = int(ells[j]) - 1;
        if (ie < weights.size()) {
          double w = weights(ie);
          double v = values(ie);
          elem(j, 0) = v;
          elem(j, 1) = w;
        } else {
          Rcpp::Rcout << "overflow: ie="<<ie<<"\n";
          break;
        }
      }
      resu.push_back(elem);
    }
  }
  return(resu);
}

//' @noRd
// //' Intersection summary for each grid cell
// //' 
// //' From the list of integer vectors which represents the intersection of grid cells and ellipses, this function returns a numeric matrix with one row per grid cell and five columns : cell gid, the number of intersecting ellipses, the weighted-averaged value of overlapping ellipses values, the sum of intersecting ellipses weights and the weighted standard deviation of overlapping ellipses values.
// //'
// //' @param cells An integer vector containing cells ids
// //' @param inter A list of integer vector containing intersections between cells and ellipses
// //' @param weights A double-precision vector containing weights of the ellipses
// //' @param values A double-precision vector containing values of the ellipses
// //' @return a numeric matrix with one row per grid cell and five columns : cell gid, the number of intersecting ellipses, the weighted-averaged value of overlapping ellipses values, the sum of intersecting ellipses weights and the weighted standard deviation of overlapping ellipses values.
// [[Rcpp::export(.parseInter_cpp)]]
Rcpp::NumericMatrix parseInter_cpp(Rcpp::NumericVector cells, Rcpp::List inter, Rcpp::DoubleVector weights, Rcpp::DoubleVector values) {
  unsigned int n = inter.size();
  Rcpp::NumericMatrix resu(n,5);
  for (unsigned int i=0; i<n; i++) {
    int gid = cells(i);
    Rcpp::List ells = inter(i);
    unsigned int ne = ells.size();
    if (ne == 0) {
      resu(i,0) = gid;
      resu(i,1) = 0;
      resu(i,1) = NA_REAL;
      resu(i,2) = NA_REAL;
      resu(i,3) = NA_REAL;
      resu(i,4) = NA_REAL;
    } else {
      // Weighted mean
      double valuesSum = 0.0;
      double weightsSum = 0.0;
      double squareSum = 0.0;
      for (unsigned int j=0; j<ne; j++) {
        int ie = int(ells[j]) - 1;
        if (ie < weights.size()) {
          double w = weights(ie);
          double v = values(ie);
          if (! (std::isnan(w) || std::isnan(v))) {
            valuesSum  += w * v;
            squareSum  += w * pow(v, 2);
            weightsSum += w;
          }
        } else {
          Rcpp::Rcout << "overflow: ie="<<ie<<"\n";
          break;
        }
      }
      double avg = valuesSum / weightsSum;
      double var = (squareSum / weightsSum) - pow(avg,2) ;
      double stdv = sqrt(var);
      resu(i,0) = gid;
      resu(i,1) = ne;
      resu(i,2) = avg;
      resu(i,3) = weightsSum;
      resu(i,4) = stdv;
    }
  }
  return(resu);
}



//' @noRd
// //' Intersection summary for each grid cell used in permutation context
// //' 
// //' From the list of integer vectors which represents the intersection of grid cells and ellipses, this function returns a numeric vector containing the weighted-averaged value of overlapping ellipses values.
// //'
// //' @param cells An integer vector containing cells ids
// //' @param inter A list of integer vector containing intersections between cells and ellipses
// //' @param weights A numeric vector containing weights of the ellipses
// //' @param values A numeric vector containing values of the ellipses
// //' @return a numeric vector with the weighted-averaged value of overlapping ellipses values.
// [[Rcpp::export(.parseInterPerm_cpp)]]
Rcpp::DoubleVector parseInterPerm_cpp(Rcpp::NumericVector cells, Rcpp::List inter, Rcpp::DoubleVector weights, Rcpp::DoubleVector values) {
  unsigned int n = inter.size();
  Rcpp::DoubleVector resu(n);
  for (unsigned int i=0; i<n; i++) {
    Rcpp::List ells = inter(i);
    unsigned int ne = ells.size();
    //Rcout << "i="<< i << "\t gid=" << gid << "\t ne=" << ne << "\n";
    if (ne == 0) {
      resu(i) = NA_REAL;
      //cout << i << "<-NA\n";
    } else {
      // Weighted mean
      double valuesSum = 0.0;
      double weightsSum = 0.0;
      for (unsigned int j=0; j<ne; j++) {
        int ie = int(ells[j]) - 1;
        if (ie < weights.size()) {
          double w = weights(ie);
          double v = values(ie);
          if (! (std::isnan(w) || std::isnan(v))) {
            valuesSum  += w * v;
            weightsSum += w;
          }
        } else {
          Rcpp::Rcout << "overflow: ie="<<ie<<"\n";
          break;
        }
      }
      resu(i) = valuesSum / weightsSum;
    }
  }
  return(resu);
}

// Function for sorting Rcomplex as if CPoint in original source
bool RcomplexSorter (Rcomplex i, Rcomplex j) { 
    //return (i.r < j.r || i.r == j.r && i.i < j.i);
    return ( (i.r < j.r) || ((i.r == j.r) && (i.i < j.i)) );
}
// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
double c_cross(const Rcomplex &O, const Rcomplex &A, const Rcomplex &B) { return (A.r - O.r) * (B.i - O.i) - (A.i - O.i) * (B.r - O.r); }

// //' @title Function convex_hull
// //' 
// //' @description Computes convex hull of cloud of points.
// //'   Implementation of Andrew's monotone chain 2D convex hull algorithm, 
// //'     modified by S. Piry for Rcomplex points, march 2018.
// //'   Asymptotic complexity: O(n log n).
// //' 
// //'   Note: the last point in the returned list is not the same as the first one (unclosed polygon).
// //' 
// //' @references https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#C++
// //' 
// //' @param P complex vector, the centroids coordinates
// //'
// //' @return a complex vector as coordinates of summits of the convex hull in counter-clockwise order.
// //'
std::vector<Rcomplex> convex_hull(Rcpp::ComplexVector P) {
  unsigned int n = P.size();
  int k = 0;
  if (n == 1) {
    std::vector<Rcomplex> H(1);
    H[0] = P[0];
    return H;
  } else {
    std::vector<Rcomplex> H(2*n);
    // Sort points lexicographically
    std::sort(P.begin(), P.end(), RcomplexSorter);
    // Build lower hull
    for (unsigned int i = 0; i < n; ++i) {
      while (k >= 2 && c_cross(H[k-2], H[k-1], P[i]) <= 0) k--;
      H[k++] = P[i];
    }
    // Build upper hull
    for (int i = n-2, t = k+1; i >= 0; i--) {
      while (k >= t && c_cross(H[k-2], H[k-1], P[i]) <= 0) k--;
      H[k++] = P[i];
    }
    H.resize(k-1);
    return H;
  }
}

// //' @title Function mkCc_cpp
// //' 
// //' @description Builds points for an circle.
// //'
// //' @param e double, the radius length
// //' @param c0 complex, coordinates of the center
// //' @param fic double vector, the N+1 angles for a quarter circle
// //'
// //' @return a complex vector as coordinates of the points
Rcpp::ComplexVector mkCc_cpp(double e, Rcomplex c0, Rcpp::DoubleVector fic) {
  int n = fic.size();
  Rcpp::ComplexVector eco(n);
  Rcpp::ComplexVector r((n-1) * 4 + 1);
  Rcomplex p, cr;
  int ip = 0;
  for (int i=0; i < n; i++) {
    p.r = e*cos(fic[i]);
    p.i = e*sin(fic[i]);
    eco[i] = p;
    cr = c0 + p;
    r[ip++] = cr;
  }
  for (int i=n-2; i > 0; i--) {
    p.r = - eco[i].r;
    p.i =   eco[i].i;
    cr = c0 + p;
    r[ip++] = cr;
  }
  for (int i=0; i < n; i++) {
    p.r = - eco[i].r;
    p.i = - eco[i].i;
    cr = c0 + p;
    r[ip++] = cr;
  }
  for (int i=n-2; i > 0; i--) {
    p.r =  eco[i].r;
    p.i = - eco[i].i;
    cr = c0 + p;
    r[ip++] = cr;
  }
  r[ip++] = r[0];
  return(r);
}

// //' @title Function mkEc_cpp
// //' 
// //' @description Builds points for an ellipse.
// //'
// //' @param a double, the half-minor axis length
// //' @param b double, the half-major axis length
// //' @param c0 complex, coordinates of the center
// //' @param rot complex, the rotation
// //' @param fic double vector, the N+1 angles for a quarter circle
// //'
// //' @return a complex vector as coordinates of the points
Rcpp::ComplexVector mkEc_cpp(double a, double b, Rcomplex c0, Rcomplex rot, Rcpp::DoubleVector fic) {
  int n = fic.size();
  Rcpp::ComplexVector eco(n);
  Rcpp::ComplexVector r((n-1) * 4 + 1);
  Rcomplex p, cr;
  int ip = 0;
  for (int i=0; i < n; ++i) {
    p.r = a*cos(fic[i]);
    p.i = b*sin(fic[i]);
    eco[i] = p;
    cr = c0 + p * rot;
    r[ip++] = cr;
  }
  for (int i=n-2; i > 0; i--) {
    p.r = - eco[i].r;
    p.i =   eco[i].i;
    cr = c0 + p * rot;
    r[ip++] = cr;
  }
  for (int i=0; i < n; ++i) {
    p.r = - eco[i].r;
    p.i = - eco[i].i;
    cr = c0 + p * rot;
    r[ip++] = cr;
  }
  for (int i=n-2; i > 0; i--) {
    p.r =   eco[i].r;
    p.i = - eco[i].i;
    cr = c0 + p * rot;
    r[ip++] = cr;
  }
  r[ip++] = r[0];
  return(r);
}

//' @noRd
// //' Ellipse-like polygon construction
// //' 
// //' Builds the convex hull of the ellipse joining two spatial points and their error circles
// //'
// //' @param r An numeric vector containing x1, y1, x2, y2, e1, e2 respectively the coordinates x,y and the error circle radius of point 1 (2).
// //' @param N The number of segments per quarter-circle
// //' @param ecc The eccentricity of the ellipse
// //' @return a numeric matrix of coordinates (x,y) in columns with one row per point of the convex hull
// //'
// [[Rcpp::export(.mkP4st_cpp)]]
Rcpp::NumericMatrix mkP4st_cpp(Rcpp::DoubleVector r, Rcpp::IntegerVector N, Rcpp::DoubleVector ecc) {
  int ic = 0;
  //     int id = r[ic++];
  double x1 = r[ic++];
  double y1 = r[ic++];
  double x2 = r[ic++];
  double y2 = r[ic++];
  double e1 = r[ic++];
  double e2 = r[ic++];
  std::complex<double> d = std::complex<double>(x1-x2, y1-y2);
  Rcomplex c0; c0.r = (x1+x2)/2.0; c0.i = (y1+y2)/2.0;
  double f = abs(d) / 2.0;
  int N0 = N[0];
  double Nf = static_cast<double>(N0);
  Rcpp::DoubleVector fic(N0+1);
  for (int i=0; i<=N0; i++) {
    fic[i] = (M_PI / 2.0L) * static_cast<double>(i) / Nf;
  }
  if (f > 0.0) { // different points
    double rot = arg(d);
    std::complex<double> crot = exp(std::complex<double>(0,1) * rot); // missing sugar...
    Rcomplex rcrot; rcrot.r=crot.real(); rcrot.i=crot.imag();
    double a = f / ecc[0];
    double b = sqrt(pow(a,2) - pow(f,2));
    Rcpp::ComplexVector ellipse = mkEc_cpp(a, b, c0, rcrot, fic);
    Rcomplex c1; c1.r=x1; c1.i=y1;
    Rcpp::ComplexVector ec1 = mkCc_cpp(e1, c1, fic);
    Rcomplex c2; c2.r=x2; c2.i=y2;
    Rcpp::ComplexVector ec2 = mkCc_cpp(e2, c2, fic);
    // Now let's build the convex hull of the three set of points
    int npE = ellipse.size();
    int npC1 = ec1.size();
    int npC2 = ec2.size();
    int ipt = 0;
    Rcpp::ComplexVector pts(npE+npC1+npC2);
    for(int i=0; i<npE;  i++) { pts[ipt++]=ellipse[i]; }
    for(int i=0; i<npC1; i++) { pts[ipt++]=ec1[i]; }
    for(int i=0; i<npC2; i++) { pts[ipt++]=ec2[i]; }
    std::vector<Rcomplex> ch1 = convex_hull(pts);
    int np = ch1.size();
    Rcpp::NumericMatrix resu(np+1, 2);
    for(int i=0; i<np; i++) {
      Rcomplex pp = ch1[i];
      resu(i,0) = pp.r;
      resu(i,1) = pp.i;
    }
    // already closed?? Not that sure...
    // close convex hull polygon
    resu(np,0) = resu(0,0);
    resu(np,1) = resu(0,1);
    return(resu);
  } else { // same point
    double e3 = (e1 > e2 ? e1 : e2);
    if (e3 > 0.0) {
      // Return error circle
      Rcomplex c1; c1.r=x1; c1.i=y1;
      Rcpp::ComplexVector ec = mkCc_cpp(e3, c1, fic);
      int np = ec.size();
      Rcpp::NumericMatrix resu(np, 2);
      for (int i=0; i < np; i++) {
        resu(i,0) = ec[i].r;
        resu(i,1) = ec[i].i;
      }
      return(resu);
    } else {
      // Return Polygon = point
      Rcpp::NumericMatrix resu(2, 2);
      resu(0,0) = r[1];
      resu(0,1) = r[2];
      return(resu);
    }
  }
}


// converts longitude,latitude to x,y,z in euclidean space on earth
std::vector<double> to_xyz_cpp(double lon1, double lat1) {
    std::vector<double> xyz(3);
    double lon = lon1 * M_PI / 180.0L;
    double lat = lat1 * M_PI / 180.0L;
    xyz[0] = EARTHRADIUS * cos(lat) * cos(lon);
    xyz[1] = EARTHRADIUS * cos(lat) * sin(lon);
    xyz[2] = EARTHRADIUS * sin(lat);
    return(xyz);
}

// converts x,y,z in euclidean space on earth to longitude,latitude
std::vector<double> from_xyz_cpp(double x, double y, double z) {
    std::vector<double> lonlat(2);
    lonlat[0] = 180.0L * atan2(y, x) / M_PI;
    lonlat[1] = 180.0L * asin(z / EARTHRADIUS) / M_PI;
    return(lonlat);
}

// computes the chord distance given a great circle distance on earth
double greatCircle2chord_cpp(double gc) {
    double d = 2.0L * EARTHRADIUS * sin(0.5L * gc / EARTHRADIUS);
    return(d);
}

// computes a great circle distance given a chord distance on earth
double chord2greatCircle_cpp(double c) {
    double alpha = asin(0.5L * c / EARTHRADIUS);
    return(2.0L * alpha * EARTHRADIUS);
}

// Summation of values with possibly huge differences
// See: https://en.wikipedia.org/wiki/Kahan_summation_algorithm

// Prevents over-optimization for next functions
#if defined(__GNUC__)
#  pragma GCC push_options
#  pragma GCC optimize ("O2")
#endif

// euclidean norm of a vector
double l2_norm_cpp(std::vector<double> const& u) {
    std::vector<double> v;
    for (unsigned int i = 0; i < u.size(); ++i)
        v.push_back(u[i]*u[i]);
    // Kahan's summation
    double accum = 0.0L;
    double c = 0.0L;
    for (unsigned int i = 0; i < v.size(); ++i) {
        volatile double y = v[i] - c;
        volatile double t = accum + y;
        volatile double z = t - accum;
        c = z - y;
        accum = t;
    }
    return(sqrt(accum));
}

// inner product of two vectors (of same length)
double in_prod_cpp(std::vector<double> const& u1, std::vector<double> const& u2) {
    if (u1.size() != u2.size())
        throw std::invalid_argument( "Vectors have different size" );
    std::vector<double> v;
    for (unsigned int i = 0; i < u1.size(); ++i)
        v.push_back(u1[i]*u2[i]);
    // Kahan's summation
    double accum = 0.0L;
    double c = 0.0L;
    for (unsigned int i = 0; i < v.size(); ++i) {
        volatile double y = v[i] - c;
        volatile double t = accum + y;
        volatile double z = t - accum;
        c = z - y;
        accum = t;
    }
    return(accum);
}

// restore default compilation options
#if defined(__GNUC__)
#  pragma GCC pop_options
#endif

// cross product of two vectors of length 3
// From: https://www.tutorialspoint.com/cplusplus-program-to-compute-cross-product-of-two-vectors
std::vector<double> cross_product_3_cpp(std::vector<double> v_A, std::vector<double> v_B) {
   std::vector<double> c_P(3);
   c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
   return(c_P);
}


//' @noRd
// //' Trilaterarion of three spheres; from original python code, see: 
// //' https://stackoverflow.com/questions/1406375/finding-intersection-points-between-3-spheres
// //' 
// //' This function is called by mkEll_s2 in MAPI_RunOnGrid.
// //' In order to build ellipse polygon, compute pairs of points as the intersection of three spheres:
// //' the earth, a sphere centered on the first focus with given radius from great circle distance,
// //' a sphere centered on the second focus with complementary radius from great circle distance
// //' such as the sum of the great circle distances is the distance between the foci divided by 
// //' eccentricity.
// //' 
//[[Rcpp::export(.trilaterate_cpp)]]
Rcpp::DoubleVector trilaterate_cpp(double p0lon, double p0lat, double p1lon, double p1lat, double gc0, double gc1) {
    Rcpp::DoubleVector resu(4);
    std::vector<double> p0_xyz = to_xyz_cpp(p0lon, p0lat);
    std::vector<double> p1_xyz = to_xyz_cpp(p1lon, p1lat);
    const std::vector<double> p3_xyz(3, 0);
    double r0 = greatCircle2chord_cpp(gc0);
    double r1 = greatCircle2chord_cpp(gc1);
    // temp1 <- P1.xyz - P0.xyz
    std::vector<double> temp1(3);
    temp1[0] = p1_xyz[0] - p0_xyz[0];
    temp1[1] = p1_xyz[1] - p0_xyz[1];
    temp1[2] = p1_xyz[2] - p0_xyz[2];
    // norm1 <- base::norm(temp1, type="2")
    double norm1 = l2_norm_cpp(temp1);
    // e_x <- temp1 / norm1
    std::vector<double> e_x(3);
    e_x[0] = temp1[0] / norm1;
    e_x[1] = temp1[1] / norm1;
    e_x[2] = temp1[2] / norm1;
    // temp2 <- P3.xyz - P0.xyz
    std::vector<double> temp2(3);
    temp2[0] = p3_xyz[0] - p0_xyz[0];
    temp2[1] = p3_xyz[1] - p0_xyz[1];
    temp2[2] = p3_xyz[2] - p0_xyz[2];
    // i <- (e_x %*% temp2)[1,1]
    double i = in_prod_cpp(e_x, temp2);
    // temp3 <- temp2 - (i * e_x)
    std::vector<double> temp3(3);
    temp3[0] = temp2[0] - i*e_x[0];
    temp3[1] = temp2[1] - i*e_x[1];
    temp3[2] = temp2[2] - i*e_x[2];
    // e_y <- temp3 / base::norm(temp3, type="2")
    double norm3 = l2_norm_cpp(temp3);
    std::vector<double> e_y(3);
    e_y[0] = temp3[0] / norm3;
    e_y[1] = temp3[1] / norm3;
    e_y[2] = temp3[2] / norm3;
    // e_z <- pracma::cross(e_x, e_y)
    std::vector<double> e_z = cross_product_3_cpp(e_x, e_y);
    // j <- (e_y %*% temp2)[1,1]
    double j = in_prod_cpp(e_y, temp2);
    // x <- (r0*r0 - r1*r1 + d*d) / (2*d)
    double x = (r0*r0 - r1*r1 + norm1*norm1) / (2.0L * norm1);
    // y <- (r0*r0 - r3*r3 - 2*i*x + i*i + j*j) / (2*j)
    double y = (r0*r0 - EARTHRADIUS*EARTHRADIUS - 2.0L*i*x + i*i + j*j) / (2.0L*j);
    // temp4 <- r0*r0 - x*x - y*y
    double temp4 = r0*r0 - x*x - y*y;
    if (isnan(temp4) || (temp4 < 0.0L)) {
        // For debug
        // Rcpp::Rcout << "The three spheres do not intersect!\n" << "p0lon=" << p0lon << " p0lat=" << p0lat << " p1lon=" << p1lon << " p1lat=" << p1lat << " gc0=" << gc0 << " gc1=" << gc1 << "\n";
        temp4 = NAN;
    }
    // z <- sqrt(temp4)
    double z = sqrt(temp4);
    // p_12_a <- P0.xyz + x*e_x + y*e_y + z*e_z
    double p12ax = p0_xyz[0] + x*e_x[0] + y*e_y[0] + z*e_z[0];
    double p12ay = p0_xyz[1] + x*e_x[1] + y*e_y[1] + z*e_z[1];
    double p12az = p0_xyz[2] + x*e_x[2] + y*e_y[2] + z*e_z[2];
    // p_12_b <- P0.xyz + x*e_x + y*e_y - z*e_z
    double p12bx = p0_xyz[0] + x*e_x[0] + y*e_y[0] - z*e_z[0];
    double p12by = p0_xyz[1] + x*e_x[1] + y*e_y[1] - z*e_z[1];
    double p12bz = p0_xyz[2] + x*e_x[2] + y*e_y[2] - z*e_z[2];
    // return(list(from.xyz(p_12_a), from.xyz(p_12_b)))
    std::vector<double> p12a = from_xyz_cpp(p12ax, p12ay, p12az);
    std::vector<double> p12b = from_xyz_cpp(p12bx, p12by, p12bz);
    resu(0) = p12a[0];
    resu(1) = p12a[1];
    resu(2) = p12b[0];
    resu(3) = p12b[1];
    return(resu);
}
