/*****************************************************************************
 * Copyright (C) 2013 x265 project
 *
 * Authors: Steve Borho <steve@borho.org>
 *          Min Chen <chenm003@163.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *
 * This program is also available under a commercial proprietary license.
 * For more information, contact us at license @ x265.com.
 *****************************************************************************/

#ifndef X265_MOTIONESTIMATE_H
#define X265_MOTIONESTIMATE_H

#include "primitives.h"
#include "reference.h"
#include "mv.h"
#include "bitcost.h"
#include "yuv.h"

namespace X265_NS {
// private x265 namespace

class MotionEstimate : public BitCost
{
protected:

    intptr_t blockOffset;					// �������׵�ַ�����֡�׵�ַ��ƫ����  
    
    int ctuAddr;
    int absPartIdx;  // part index of PU, including CU offset within CTU

    int searchMethod;					// ME�����㷨  
    int subpelRefine;					// submeǿ��  

	int blockwidth;						// ��ǰ������Ŀ��
	int blockheight;					// ��x265�������ʱû���κ���;  

    pixelcmp_t sad;						// ���ڼ�����SADֵ (src, srcStride, dst, dstStride,);  
    pixelcmp_x3_t sad_x3;				// ͬʱ����3��MV��Ӧ��3��SADֵ			template<int lx, int ly> void sad_x3(const pixel* pix1, const pixel* pix2, const pixel* pix3, const pixel* pix4, intptr_t frefstride, int32_t* res)  
    pixelcmp_x4_t sad_x4;				// ͬʱ����4��MV��Ӧ��4��SADֵ template<int lx, int ly> oid sad_x4(const pixel* pix1, const pixel* pix2, const pixel* pix3, const pixel* pix4, const pixel* pix5, intptr_t frefstride, int32_t* res)  
    pixelcmp_ads_t ads;					// ����SATDֵ��������̿��Բ鿴satd_4x4����(fencIntra, cuSize, prediction, cuSize)
    pixelcmp_t satd;
    pixelcmp_t chromaSatd;

    MotionEstimate& operator =(const MotionEstimate&);

public:

    static const int COST_MAX = 1 << 28;

    uint32_t* integral[INTEGRAL_PLANE_NUM];
    Yuv fencPUYuv;						// ��������Ļ���,��СΪ64x64���������������copy���˻���  
    int partEnum;
    bool bChromaSATD;					// �Ƿ����chroma������satd��ֻ����subpelRefine����2ʱ���ڷ�����MEʱ�Ż����chroma��satd  

	/** ��������             ����ʼ��ME��searchMethod Ĭ��hex��subme Ĭ��2
	* \����                  ��null * */
    MotionEstimate();
    ~MotionEstimate();

    static void initScales();
    static int hpelIterationCount(int subme);

	/** ��������   �� ��ʼ�������㷨��������������Ļ���
	/*\����  method�� ��������
	/*\����  refine�� submeǿ��
	* \����     csp:  ͼ���ʽ
	* \����        �� null */
    void init(int csp);

    /* Methods called at slice setup */
	/* Methods called at slice setup */
	/** ��������   �� ����me��Ӧ��asm������copy�����������ݵ���������Ļ���
	/*\����   fencY�� ��ǰ����֡��֡�׵�ַ
	/*\����  stride�� ԭʼ֡����
	/*\����  offset�� ��ǰ�������׵�ַ�����֡�׵�ַ��ƫ����
	/*\����  pwidth�� ��ǰ������Ŀ��
	/*\���� pheight�� ��ǰ������ĸ߶�
	* \����        �� null */
    void setSourcePU(pixel *fencY, intptr_t stride, intptr_t offset, int pwidth, int pheight, const int searchMethod, const int subpelRefine);
    void setSourcePU(const Yuv& srcFencYuv, int ctuAddr, int cuPartIdx, int puPartIdx, int pwidth, int pheight, const int searchMethod, const int subpelRefine, bool bChroma);

    /* buf*() and motionEstimate() methods all use cached fenc pixels and thus
     * require setSourcePU() to be called prior. */
	/** ��������   �� ���㵱ǰ����ο����SATDֵ
	/*\����    fref�� �ο����׵�ַ
	/*\����  stride�� �ο��鲽������
	* \����        �� ��ǰ����ο����SATDֵ */
    inline int bufSAD(const pixel* fref, intptr_t stride)  { return sad(fencPUYuv.m_buf[0], FENC_STRIDE, fref, stride); }

    inline int bufSATD(const pixel* fref, intptr_t stride) { return satd(fencPUYuv.m_buf[0], FENC_STRIDE, fref, stride); }

    inline int bufChromaSATD(const Yuv& refYuv, int puPartIdx)
    {
        return chromaSatd(refYuv.getCbAddr(puPartIdx), refYuv.m_csize, fencPUYuv.m_buf[1], fencPUYuv.m_csize) +
               chromaSatd(refYuv.getCrAddr(puPartIdx), refYuv.m_csize, fencPUYuv.m_buf[2], fencPUYuv.m_csize);
    }

    int motionEstimate(ReferencePlanes* ref, const MV & mvmin, const MV & mvmax, const MV & qmvp, int numCandidates, const MV * mvc, int merange, MV & outQMv, pixel *srcReferencePlane = 0);

    int subpelCompare(ReferencePlanes* ref, const MV &qmv, pixelcmp_t);

protected:
	/**
	* �������ܣ�����ME����
	/*  ���÷�Χ        ֻ��MotionEstimate::motionEstimate�����б�����
	* \����  ref        �ο�֡
	* \����  mvmin      �����ʵ��������Χ����߽���ϱ߽磩
	* \����  mvmax      �����ʵ��������Χ���±߽���ұ߽磩
	* \����  bmv        ��AMVP�õ���Ԥ��MV�����������ŵ�MV
	* \����  bcost      Ԥ��MV��Ӧ��cost�����������ŵ�cost
	* \����  bPointNr   �������ŵ�MV��Ӧ��λ�ñ�ţ���λ�ñ��������ME������ģ���б��
	* \����  bDistance  �������ŵ�MV��Ӧ�Ĳ���
	* \����  earlyExitIters ��ǰ�����ĵ�������
	* \����  merange    �����ME������Χ
	*/
    inline void StarPatternSearch(ReferencePlanes *ref,
                                  const MV &       mvmin,
                                  const MV &       mvmax,
                                  MV &             bmv,
                                  int &            bcost,
                                  int &            bPointNr,
                                  int &            bDistance,
                                  int              earlyExitIters,
                                  int              merange);
};
}

#endif // ifndef X265_MOTIONESTIMATE_H
