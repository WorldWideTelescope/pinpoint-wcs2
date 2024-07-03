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

#include <iostream>
#include <fstream>
#include <string>

#include <QDebug>
#include <QImage>
#include <QFile>

#include <fitsio.h>

#include "version.h"
#include "ExportWCS.h"

// Ensure XMP templates are instantiated
#include "XMP.incl_cpp"

// Define the AVM namespace
#define kXMP_NS_AVM "http://www.communicatingastronomy.org/avm/1.0/"
#define kXMP_NS_CXC "http://pinpoint.cfa.harvard.edu"


ExportWCS::ExportWCS(QString *f, QPixmap *p, ComputeWCS *cwcs)
{
    filename = f;
    pixmap = p;
    computewcs = cwcs;
    fitsexport = false;
}


ExportWCS::~ExportWCS()
{}


void
ExportWCS::setWCS(struct WorldCoor *w)
{
    wcs = w;
}


void
ExportWCS::clearWCS()
{
    wcs = NULL;
}


void
ExportWCS::exportFITS()
{
    // Initialize variables
    fitsfile *fptr;
    int status = 0;
    int ii, jj;
    long fpixel, nelements;
    unsigned short *imagedata[pixmap->height()];
    QImage im = pixmap->toImage();

    // Prompt user for filename
    saveas = QFileDialog::getSaveFileName(NULL, "Export FITS Image", *filename+"_ppwcs.fits", tr("Images(*.fit *.fits)"));

    // Check that a filename is chosen
    if (saveas.isEmpty())
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }

    // Initialize FITS image parameters
    int bitpix = BYTE_IMG;
    long naxis = 2;
    long naxes[2] = {pixmap->width(), pixmap->height()};

    // Allocate memory for the entire image
    imagedata[0] = (unsigned short *) malloc(naxes[0] * naxes[1] * sizeof(unsigned short));
    if (imagedata[0] == NULL)
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }

    // Initialize pointers to the start of each row of the image
    for(ii=1; ii<naxes[1]; ii++)
        imagedata[ii] = imagedata[ii-1] + naxes[0];

    // Remove file is it already exists
    remove(saveas.toStdString().c_str());

    // Attempt to create the FITS file (writes to disk)
    if (fits_create_file(&fptr, saveas.toStdString().c_str(), &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }

    // Attempt to create the image HDU
    if (fits_create_img(fptr,  bitpix, naxis, naxes, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }

    // Initialize the values in the image using the QPixmaps
    for (jj = 0; jj < naxes[1]; jj++)
        for (ii = 0; ii < naxes[0]; ii++)
            imagedata[jj][ii] = qGray(im.pixel(ii, naxes[1]-jj-1));

    // First pixel and number of pixels to write
    fpixel = 1;
    nelements = naxes[0] * naxes[1];

    // Write the array to the FITS image
    if (fits_write_img(fptr, TUSHORT, fpixel, nelements, imagedata[0], &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }

    // Free memory
    free(imagedata[0]);

    // Initialize variables for the FITS header
    char xtension[] = "IMAGE";
    char origin[] = "PinpointWCS by the Chandra X-ray Center";
    char wcsname[] = "Primary WCS";
    char ctype1[] = "RA---TAN";
    char ctype2[] = "DEC--TAN";
    char cunit[] = "deg";

    // TSTRING, TLOGICAL (== int), TBYTE, TSHORT, TUSHORT, TINT, TUINT, TLONG, TLONGLONG, TULONG, TFLOAT, TDOUBLE
    if (fits_update_key(fptr, TSTRING, "XTENSION", &xtension, NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "ORIGIN", &origin, NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TINT, "WCSAXES", &naxis, NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "WCSNAME", &wcsname, NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }

    // Format values to specific digits
    QString equinox = QString("%1").arg(wcs->equinox, 0, 'f', 1);
    QString crpix1 = QString("%1").arg(wcs->xrefpix, 0, 'f', 11);
    QString crval1 = QString("%1").arg(wcs->xref, 0, 'f', 11);
    QString crpix2 = QString("%1").arg(wcs->yrefpix, 0, 'f', 11);
    QString crval2 = QString("%1").arg(wcs->yref, 0, 'f', 11);
    QString cd11 = QString("%1").arg(wcs->cd[0], 0, 'f', 11);
    QString cd12 = QString("%1").arg(wcs->cd[1], 0, 'f', 11);
    QString cd21 = QString("%1").arg(wcs->cd[2], 0, 'f', 11);
    QString cd22 = QString("%1").arg(wcs->cd[3], 0, 'f', 11);

    if (fits_update_key(fptr, TSTRING, "EQUINOX", (void*) equinox.toStdString().c_str(), NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "RADESYS", &(wcs->radecsys), NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CTYPE1", &ctype1, NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CRPIX1", (void*) crpix1.toStdString().c_str(), NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CRVAL1", (void*) crval1.toStdString().c_str(), NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CUNIT1", &cunit, NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CTYPE2", &ctype2, NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CRPIX2", (void*) crpix2.toStdString().c_str(), NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CRVAL2", (void*) crval2.toStdString().c_str(), NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CUNIT2", &cunit, NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CD1_1", (void*) cd11.toStdString().c_str(), NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CD1_2", (void*) cd21.toStdString().c_str(), NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CD2_1", (void*) cd12.toStdString().c_str(), NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }
    if (fits_update_key(fptr, TSTRING, "CD2_2", (void*) cd22.toStdString().c_str(), NULL, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }

    // Write comments
    if (fits_write_comment(fptr, "World Coordinate System computed using PinpointWCS by the Chandra X-ray Center.  PinpointWCS is developed and maintained by Amit Kapadia (CfA) akapadia@cfa.harvard.edu.", &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }

    // Close FITS file
    if (fits_close_file(fptr, &status))
    {
        fitsexport = false;
        emit exportResults(fitsexport);
        return;
    }

    fitsexport = true;
    emit exportResults(true);
}


void
ExportWCS::exportAVMClean()
{
    exportAVM(false);
}


void
ExportWCS::exportAVMDetailed()
{
    exportAVM(true);
}


void
ExportWCS::exportXMP()
{
    SXMPMeta avm;
    bool success = false;

    if (!constructAVM(&avm, false, NULL)) {
        saveas = QFileDialog::getSaveFileName(
            NULL,
            "Export XMP Packet",
            *filename + ".xmp",
            tr("Metadata (*.xmp))")
        );

        if (!saveas.isEmpty())
        {
            std::string xmpstr;
            avm.SerializeToBuffer(&xmpstr);

            std::ofstream outFile;
            outFile.open(saveas.toStdString().c_str());
            outFile << xmpstr;
            outFile.close();
            success = true;
        }

        SXMPFiles::Terminate();
        SXMPMeta::Terminate();
    }

    emit exportResults(success);
}


void
ExportWCS::exportAVM(bool detailed)
{
    SXMPMeta avm;
    SXMPFiles epoimage;
    bool success = false;

    if (!constructAVM(&avm, detailed, &epoimage)) {
        if (epoimage.CanPutXMP(avm)) {
            epoimage.PutXMP(avm);
            success = true;
        } else
            success = false;

        epoimage.CloseFile();
        SXMPFiles::Terminate();
        SXMPMeta::Terminate();
    }

    emit exportResults(success);
}


// Returns false on success, true on error.
//
// If this function succeeds, the caller must call SXMPFiles::Terminate() and
// SXMPMeta::Terminate() after finishing using the AVM data.
//
// If `out_epoimage` is NULL, the opened XMP handle to the EPO image will be
// closed when this function finishes. Otherwise, the file will *not* be closed,
// and the file handle will be stored in the pointer location. The caller must
// close the handle before calling the Terminate() functions.
bool
ExportWCS::constructAVM(SXMPMeta *avm, bool detailed, SXMPFiles *out_epoimage)
{
    bool error = true;
    std::string f = filename->toStdString();

    if (!SXMPMeta::Initialize())
        return error;

    XMP_OptionBits init_options = 0;
#if UNIX_ENV
    init_options |= kXMPFiles_ServerMode;
#endif

    if (!SXMPFiles::Initialize(init_options)) {
        SXMPMeta::Terminate();
        return error;
    }

    try
    {
        XMP_OptionBits open_options = kXMPFiles_OpenForUpdate | kXMPFiles_OpenUseSmartHandler;
        bool ok;
        SXMPFiles epoimage;

        ok = epoimage.OpenFile(f, kXMP_UnknownFile, open_options);
        if (!ok)
        {
            open_options = kXMPFiles_OpenForUpdate | kXMPFiles_OpenUsePacketScanning;
            ok = epoimage.OpenFile(f, kXMP_UnknownFile, open_options);
        }

        if (ok) {
            epoimage.GetXMP(avm);

            std::string avmprefix;
            SXMPMeta::RegisterNamespace(kXMP_NS_AVM, "avm", &avmprefix);

            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.CoordinateFrame");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.Equinox");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.ReferenceValue");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.ReferenceDimension");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.ReferencePixel");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.Scale");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.Rotation");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.CoordsystemProjection");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.Quality");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.Notes");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.FITSheader");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Spatial.CDMatrix");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Publisher.MetadataDate");
            avm->DeleteProperty(kXMP_NS_AVM, "avm:Publisher.MetadataVersion");

            QString equinox = QString("%1").arg(wcs->equinox, 0, 'f', 1);
            QString crval1 = QString("%1").arg(computewcs->crval(0), 0, 'f', 11);
            QString crval2 = QString("%1").arg(computewcs->crval(1), 0, 'f', 11);
            QString crpix1 = QString("%1").arg(computewcs->crpix(0), 0, 'f', 11);
            QString crpix2 = QString("%1").arg(computewcs->epo_height - computewcs->crpix(1), 0, 'f', 11);
            QString scale1 = QString("%1").arg(-1*computewcs->scale, 0, 'f', 11);
            QString scale2 = QString("%1").arg(computewcs->scale, 0, 'f', 11);
            QString orientation = QString("%1").arg(computewcs->orientation, 0, 'f', 11);
            QString width = QString("%1").arg(computewcs->epo_width);
            QString height = QString("%1").arg(computewcs->epo_height);
            QString spatialnotes = QString("World Coordinate System resolved using PinpointWCS %1 revision %2 by the Chandra X-ray Center").arg(VERSION).arg(REVISION);

            if (detailed)
            {
                QString data = QString("\n\n%1\t\t%2\t\t%3\t\t%4\n").arg("FITS X").arg("FITS Y").arg("EPO X").arg("EPO Y");
                spatialnotes.append(data);

                for (int i = 0; i < computewcs->refCoords->size(); i++)
                {
                    QString data = QString("%1\t\t%2\t\t%3\t\t%4\n").arg(computewcs->refCoords->at(i).x(), 0, 'f', 2).arg(computewcs->refCoords->at(i).y(), 0, 'f', 2).arg(computewcs->epoCoords->at(i).x(), 0, 'f', 2).arg(computewcs->epoCoords->at(i).y(), 0, 'f', 2);
                    spatialnotes.append(data);
                }

                // Get the center pixel (for STScI)
                QString center_x = QString("\n%1").arg(computewcs->epo_width/2., 0, 'f', 2);
                QString center_y = QString("%1").arg(computewcs->epo_height/2., 0, 'f', 2);
                QString center_ra = QString("FIXME"); // QString("%1").arg(computewcs->centerRA, 0, 'f', 11);
                QString center_dec = QString("FIXME"); // QString("%1").arg(computewcs->centerDec, 0, 'f', 11);

                QString centerpix = QString("\nCenter Pixel Coordinates:%1\t%2\n%3\t%4").arg(center_x).arg(center_ra).arg(center_y).arg(center_dec);
                spatialnotes.append(centerpix);
            }

            XMP_OptionBits itemOptions;
            itemOptions = kXMP_PropValueIsArray | kXMP_PropArrayIsOrdered;

            avm->SetProperty(kXMP_NS_AVM, "avm:Spatial.CoordinateFrame", "ICRS", 0);
            avm->SetProperty(kXMP_NS_AVM, "avm:Spatial.Equinox", equinox.toStdString(), 0);
            avm->AppendArrayItem(kXMP_NS_AVM, "avm:Spatial.ReferenceValue", itemOptions, crval1.toStdString());
            avm->AppendArrayItem(kXMP_NS_AVM, "avm:Spatial.ReferenceValue", itemOptions, crval2.toStdString());
            avm->AppendArrayItem(kXMP_NS_AVM, "avm:Spatial.ReferenceDimension", itemOptions, width.toStdString());
            avm->AppendArrayItem(kXMP_NS_AVM, "avm:Spatial.ReferenceDimension", itemOptions, height.toStdString());
            avm->AppendArrayItem(kXMP_NS_AVM, "avm:Spatial.ReferencePixel", itemOptions, crpix1.toStdString());
            avm->AppendArrayItem(kXMP_NS_AVM, "avm:Spatial.ReferencePixel", itemOptions, crpix2.toStdString());
            avm->AppendArrayItem(kXMP_NS_AVM, "avm:Spatial.Scale", itemOptions, scale1.toStdString());
            avm->AppendArrayItem(kXMP_NS_AVM, "avm:Spatial.Scale", itemOptions, scale2.toStdString());
            avm->SetProperty(kXMP_NS_AVM, "avm:Spatial.Rotation", orientation.toStdString(), 0);
            avm->SetProperty(kXMP_NS_AVM, "avm:Spatial.CoordsystemProjection", "TAN", 0);
            avm->SetProperty(kXMP_NS_AVM, "avm:Spatial.Quality", "Full", 0);
            avm->SetLocalizedText(kXMP_NS_AVM, "avm:Spatial.Notes", "x-default", "x-default", spatialnotes.toStdString(), 0);
            // avm->SetProperty(kXMP_NS_AVM, "avm:Spatial.FITSheader", ...", 0);

            XMP_DateTime updatedTime;
            SXMPUtils::CurrentDateTime(&updatedTime);
            avm->SetProperty_Date(kXMP_NS_AVM, "avm:MetadataDate", updatedTime, 0);
            avm->SetProperty(kXMP_NS_AVM, "avm:MetadataVersion", AVM_VERSION, 0);

            if (out_epoimage == nullptr)
                epoimage.CloseFile();
            else
                *out_epoimage = epoimage;

            error = false;
        }
    }
    catch (XMP_Error &e)
    {
        std::cerr << "XMP Error: " << e.GetErrMsg() << std::endl;
    }

    return error;
}
