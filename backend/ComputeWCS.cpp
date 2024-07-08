/*
 *  PinpointWCS is developed by the Chandra X-ray Center
 *  Education and Public Outreach Group
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

using namespace std;

#include <cmath>
#include <iostream>

#include <Eigen/LU>
#include <QDebug>

#include "ComputeWCS.h"
#include "PinpointWCSUtils.h"


ComputeWCS::ComputeWCS(
    QList<QPointF> *ref,
    QList<QPointF> *epo,
    struct WorldCoor *refWCS,
    int epo_width,
    int epo_height
) {
    // Input coordinate convention:
    // ref (FITS): see fitsPixelCoordinate() below
    // epo (JPEG): JPEG-like, Y = 0 is the top.

    refCoords = ref;
    epoCoords = epo;
    referenceWCS = refWCS;
    epoWCS = false;
    this->epo_width = epo_width;
    this->epo_height = epo_height;
    downsample_factor = 1;
    rms_x = 0.0;
    rms_y = 0.0;
}


ComputeWCS::~ComputeWCS() {}


QPointF
ComputeWCS::fitsPixelCoordinate(int i)
{
    // Convert our "refCoords" to true FITS pixel coordinates. In FITS, the
    // center of the first pixel of a 2D image has coordinate (1,1); the pixel's
    // edges go from 0.5 to 1.5. Meanwhile, we're passed coordinates where the
    // pixel centers are half-integers, *and* there has been a flip of the Y
    // axis to match the usual JPEG-like ordering: y = 0.5 denotes the top row
    // of pixels in the image, which is "really" the row with `y = image_height`
    // according to FITS.
    //
    // Furthermore, if downsampling is active, the coordinates that we have been
    // given are all referenced to a version of the image that has been
    // downsampled by an integer factor. This should be straightforward to
    // handle, but there appears to be something wonky about how the coordinates
    // are computed elsewhere in the program. (For instance, the circle
    // indicating the selected FITS point is not precisely centered on the
    // selected point!) In order to get accurate results for different
    // downsample factors, I need to use the following expressions.
    QPointF raw = refCoords->at(i);
    float x = downsample_factor * (raw.x() - 0.5) + 1;
    float y = downsample_factor * (raw.y() - 0.5) + 1;
    return QPointF(x, referenceWCS->nypix - y);
}

void
ComputeWCS::computeTargetWCS()
{
    int numPoints = refCoords->size();

    if (numPoints > 0 && epoCoords->last() == QPointF(-1, -1)) {
        numPoints--;
    }

    if (numPoints < 3) {
        epoWCS = false;
        emit nowcs();
        return;
    }

    // OK, we can solve! First step: compute the affine transformation matrix.
    //
    // TODO: ensure that pixel-numbering conventions are all correct.

    Vector3d basis = Vector3d::Zero();
    Matrix3d matrix = Matrix3d::Zero();
    Vector3d xvector = Vector3d::Zero();
    Vector3d yvector = Vector3d::Zero();

    for (int i = 0; i < numPoints; i++)
    {
        QPointF p_fits = fitsPixelCoordinate(i);
        QPointF p_epo = epoCoords->at(i);

        std::cout << "P: " << p_fits.x() << "\t" << p_fits.y() << "\t" << p_epo.x() << "\t" << p_epo.y() << "\n";

        basis << p_epo.x(), p_epo.y(), 1;
        matrix += basis * basis.transpose();
        xvector += p_fits.x() * basis;
        yvector += p_fits.y() * basis;
    }

    Vector3d xcoeff = matrix.lu().solve(xvector);
    Vector3d ycoeff = matrix.lu().solve(yvector);
    Matrix3d epo2fits = Matrix3d::Identity();
    epo2fits.row(0) = xcoeff;
    epo2fits.row(1) = ycoeff;

    Matrix3d fits2epo = epo2fits.inverse().eval();

    // Construct information for the desired WCS.
    //
    // CRVAL is the same as the reference image.

    crval(0) = referenceWCS->crval[0];
    crval(1) = referenceWCS->crval[1];

    // CRPIX is the pixel from the reference image, mapped
    // to the EPO image.

    Vector3d vwork = Vector3d();
    vwork << referenceWCS->crpix[0], referenceWCS->crpix[1], 1;
    vwork = fits2epo * vwork;
    crpix(0) = vwork[0];
    crpix(1) = vwork[1];

    // CD matrix concatenates the transforms. Since the CD matrix is of
    // contravariant type, we have to multiply it by epo2fits, not fits2epo.

    Matrix2d mwork = Matrix2d();
    cdmatrix << referenceWCS->cd[0], referenceWCS->cd[1], referenceWCS->cd[2], referenceWCS->cd[3];
    mwork << epo2fits(0,0), epo2fits(0,1), epo2fits(1,0), epo2fits(1,1);
    cdmatrix = mwork * cdmatrix;

    std::cout << "crval: " << crval << std::endl;
    std::cout << "crpix: " << crpix << std::endl;
    std::cout << "CD: " << cdmatrix << std::endl;

    // Done.

    epoWCS = true;
    emit wcs();
}


struct WorldCoor *
ComputeWCS::makeTargetWCS()
{
    struct WorldCoor *targetWCS = wcskinit(
        epo_width,
        epo_height,
        (char *) "RA---TAN",
        (char *) "DEC--TAN",
        crpix(0),
        crpix(1),
        crval(0),
        crval(1),
        cdmatrix.data(),
        0.0, // cdelt1, unused if CD matrix is defined
        0.0, // cdelt2, ditto
        0.0, // crota, ditto
        referenceWCS->equinox,
        referenceWCS->epoch
    );

    wcsoutinit(targetWCS, (char *) "FK5");
    return targetWCS;
}


void ComputeWCS::computeResiduals(int numPoints)
{
#if 0
    // Initialize some variables
    double sumn = 0;
    double sumx2 = 0;
    double sumy2 = 0;

    // Loop over the pairs stored in the data model
    for (int ii=0; ii < numPoints; ii++)
    {
        // Store the coordinate pairs
        QPointF point1 = refCoords->at(ii);
        QPointF point2 = epoCoords->at(ii);

        // Map EPO coordinates to FITS coordinates
        Vector2d fit = epoToFits(&point2);

        // Compute residuals
        sumn += 1;
        sumx2 += pow(point1.x() - fit[0], 2);
        sumy2 += pow(point1.y() - fit[1], 2);
    }

    rms_x = sqrt(sumx2 / sumn);
    rms_y = sqrt(sumy2 / sumn);
#else
    rms_x = rms_y = 0;
#endif
}


Vector2d
ComputeWCS::fitsToEpo(double x, double y)
{
#if 0
    // TODO: Check that bd-af is nonzero
    Vector2d epoCoord;
    float x0, y0;

    // Create matrix of coeffients
    Matrix2d m;
    m << xcoeff(0), xcoeff(1), ycoeff(0), ycoeff(1);

    // Take the inverse
    m = m.inverse().eval();
    Vector3d xinverse;
    Vector3d yinverse;
    xinverse << m(0,0), m(0,1), xcoeff[2];
    yinverse << m(1,0), m(1,1), ycoeff[2];

    // Find the difference
    x0 = x - xinverse[2];
    y0 = y - yinverse[2];

    epoCoord << xinverse[0] * x0 + xinverse[1] * y0, yinverse[0] * x0 + yinverse[1] * y0;
    return epoCoord;
#else
    return Vector2d::Zero();
#endif
}


QPointF
ComputeWCS::fitsToEpo(QPointF *p_fits)
{
    Vector2d p_epo = fitsToEpo(p_fits->x(), p_fits->y());
    return QPointF(p_epo[0], p_epo[1]);
}


void
ComputeWCS::setDownsampleFactor(int factor)
{
    downsample_factor = factor;
}
