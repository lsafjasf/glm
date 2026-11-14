#define GLM_ENABLE_EXPERIMENTAL

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtx/matrix_query.hpp>
#include <glm/gtx/quaternion.hpp>

#include <cmath>
#include <cstdio>
#include <limits>

namespace precision_regression
{
	template<typename T>
	using Mat3 = glm::mat<3, 3, T, glm::defaultp>;

	template<typename T>
	using Mat4 = glm::mat<4, 4, T, glm::defaultp>;

	template<typename T>
	using Vec3 = glm::vec<3, T, glm::defaultp>;

	template<typename T>
	using Vec4 = glm::vec<4, T, glm::defaultp>;

	template<typename T>
	using Quat = glm::qua<T, glm::defaultp>;

	template<typename T>
	struct Thresholds;

	// See test/precision/README.md for the threshold table and fixed inputs.
	template<>
	struct Thresholds<float>
	{
		static constexpr double ProjectionAbs = 2.0e-7;
		static constexpr double ProjectionRel = 2.0e-6;
		static constexpr double ExtremeProjectionAbs = 1.0e-10;
		static constexpr double ExtremeProjectionRel = 2.0e-5;
		static constexpr double Orthogonality = 5.0e-6;
		static constexpr double Determinant = 5.0e-6;
		static constexpr double Decomposition = 3.0e-6;
		static constexpr double TinyScaleAbs = 1.0e-12;
		static constexpr double TinyScaleRel = 3.0e-6;
		static constexpr double EulerAngle = 2.0e-6;
		static constexpr double Quaternion = 6.0e-6;
		static constexpr double QuaternionAngle = 6.0e-6;
	};

	template<>
	struct Thresholds<double>
	{
		static constexpr double ProjectionAbs = 2.0e-14;
		static constexpr double ProjectionRel = 2.0e-14;
		static constexpr double ExtremeProjectionAbs = 1.0e-18;
		static constexpr double ExtremeProjectionRel = 2.0e-14;
		static constexpr double Orthogonality = 2.0e-14;
		static constexpr double Determinant = 2.0e-14;
		static constexpr double Decomposition = 2.0e-14;
		static constexpr double TinyScaleAbs = 1.0e-20;
		static constexpr double TinyScaleRel = 2.0e-14;
		static constexpr double EulerAngle = 2.0e-14;
		static constexpr double Quaternion = 3.0e-14;
		static constexpr double QuaternionAngle = 3.0e-14;
	};

	static int Errors = 0;

	static void checkMetric(
		char const* test,
		char const* metric,
		double measured,
		double limit)
	{
		if(measured <= limit)
			return;

		std::printf(
			"FAIL %s: %s measured=%.17g limit=%.17g\n",
			test, metric, measured, limit);
		++Errors;
	}

	static void checkTrue(char const* test, char const* condition, bool result)
	{
		if(result)
			return;

		std::printf("FAIL %s: %s\n", test, condition);
		++Errors;
	}

	template<typename T>
	static bool hasNonFinite(Mat4<T> const& matrix)
	{
		for(glm::length_t column = 0; column < 4; ++column)
		for(glm::length_t row = 0; row < 4; ++row)
			if(!std::isfinite(matrix[column][row]))
				return true;
		return false;
	}

	template<typename T>
	static bool hasNonFinite(Quat<T> const& quaternion)
	{
		return
			!std::isfinite(quaternion.w) ||
			!std::isfinite(quaternion.x) ||
			!std::isfinite(quaternion.y) ||
			!std::isfinite(quaternion.z);
	}

	template<typename T>
	static double matrixError(
		Mat4<T> const& actual,
		Mat4<T> const& expected,
		double absoluteTolerance,
		double relativeTolerance)
	{
		double result = 0.0;
		for(glm::length_t column = 0; column < 4; ++column)
		for(glm::length_t row = 0; row < 4; ++row)
		{
			double const actualValue = static_cast<double>(actual[column][row]);
			double const expectedValue = static_cast<double>(expected[column][row]);
			double const denominator =
				absoluteTolerance + relativeTolerance * std::abs(expectedValue);
			double const error = std::abs(actualValue - expectedValue) / denominator;
			if(!std::isfinite(error))
				return std::numeric_limits<double>::infinity();
			result = glm::max(result, error);
		}
		return result;
	}

	template<typename T>
	static double quaternionError(
		Quat<T> const& actual,
		Quat<T> const& expected,
		double absoluteTolerance,
		double relativeTolerance)
	{
		double result = 0.0;
		T const values[] = {actual.w, actual.x, actual.y, actual.z};
		T const references[] = {expected.w, expected.x, expected.y, expected.z};
		for(int i = 0; i < 4; ++i)
		{
			double const denominator =
				absoluteTolerance + relativeTolerance * std::abs(static_cast<double>(references[i]));
			double const error = std::abs(
				static_cast<double>(values[i]) - static_cast<double>(references[i])) / denominator;
			if(!std::isfinite(error))
				return std::numeric_limits<double>::infinity();
			result = glm::max(result, error);
		}
		return result;
	}

	static double scalarError(
		double actual,
		double expected,
		double absoluteTolerance,
		double relativeTolerance)
	{
		double const denominator =
			absoluteTolerance + relativeTolerance * std::abs(expected);
		double const error = std::abs(actual - expected) / denominator;
		return std::isfinite(error) ? error : std::numeric_limits<double>::infinity();
	}

	template<typename T>
	static T rotationAngle(Quat<T> const& quaternion)
	{
		T const vectorLength =
			glm::length(Vec3<T>(quaternion.x, quaternion.y, quaternion.z));
		return static_cast<T>(2) * static_cast<T>(std::atan2(
			static_cast<double>(vectorLength),
			static_cast<double>(glm::abs(quaternion.w))));
	}

	template<typename T>
	static void printThresholds(char const* name)
	{
		std::printf(
			"%s thresholds: projection=%.3g abs/%.3g rel, "
			"extreme=%.3g abs/%.3g rel, orthogonality=%.3g, determinant=%.3g, "
			"decompose=%.3g, euler=%.3g rad, quaternion=%.3g\n",
			name,
			Thresholds<T>::ProjectionAbs, Thresholds<T>::ProjectionRel,
			Thresholds<T>::ExtremeProjectionAbs, Thresholds<T>::ExtremeProjectionRel,
			Thresholds<T>::Orthogonality, Thresholds<T>::Determinant,
			Thresholds<T>::Decomposition, Thresholds<T>::EulerAngle,
			Thresholds<T>::Quaternion);
	}

	template<typename T>
	static void testViewMatrices()
	{
		char const* test = "view-hand";
		Vec3<T> const eye(static_cast<T>(1.3), static_cast<T>(-2.2), static_cast<T>(4.7));
		Vec3<T> const center(static_cast<T>(0.4), static_cast<T>(1.1), static_cast<T>(0.2));
		Vec3<T> const up(static_cast<T>(0.0), static_cast<T>(1.0), static_cast<T>(0.0));
		T const distance = glm::length(center - eye);

		Mat4<T> const right = glm::lookAtRH(eye, center, up);
		Mat4<T> const left = glm::lookAtLH(eye, center, up);
		Vec4<T> const rightCenter = right * Vec4<T>(center, static_cast<T>(1));
		Vec4<T> const leftCenter = left * Vec4<T>(center, static_cast<T>(1));

		double const positionTolerance =
			Thresholds<T>::Orthogonality * (1.0 + static_cast<double>(distance));
		checkMetric(test, "RH center x", std::abs(static_cast<double>(rightCenter.x)), positionTolerance);
		checkMetric(test, "RH center y", std::abs(static_cast<double>(rightCenter.y)), positionTolerance);
		checkMetric(
			test, "RH center z",
			std::abs(static_cast<double>(rightCenter.z) + static_cast<double>(distance)),
			positionTolerance);
		checkMetric(test, "LH center x", std::abs(static_cast<double>(leftCenter.x)), positionTolerance);
		checkMetric(test, "LH center y", std::abs(static_cast<double>(leftCenter.y)), positionTolerance);
		checkMetric(
			test, "LH center z",
			std::abs(static_cast<double>(leftCenter.z) - static_cast<double>(distance)),
			positionTolerance);

		Mat4<T> const views[] = {right, left};
		for(int viewIndex = 0; viewIndex < 2; ++viewIndex)
		{
			Mat3<T> linear;
			for(glm::length_t column = 0; column < 3; ++column)
			for(glm::length_t row = 0; row < 3; ++row)
				linear[column][row] = views[viewIndex][column][row];

			for(glm::length_t columnA = 0; columnA < 3; ++columnA)
			for(glm::length_t columnB = 0; columnB < 3; ++columnB)
			{
				double const expected = columnA == columnB ? 1.0 : 0.0;
				double const measured = static_cast<double>(
					glm::dot(linear[columnA], linear[columnB]));
				checkMetric(
					viewIndex == 0 ? "view-rh" : "view-lh",
					"orthonormal basis dot",
					std::abs(measured - expected),
					Thresholds<T>::Orthogonality);
			}

			checkMetric(
				viewIndex == 0 ? "view-rh" : "view-lh",
				"linear determinant",
				std::abs(static_cast<double>(glm::determinant(linear)) - 1.0),
				Thresholds<T>::Determinant);
			checkMetric(
				viewIndex == 0 ? "view-rh" : "view-lh",
				"affine determinant",
				std::abs(static_cast<double>(glm::determinant(views[viewIndex])) - 1.0),
				Thresholds<T>::Determinant);
		}
	}

	template<typename T>
	static Mat4<T> makePerspectiveReference(T depth22, T depth23, T depth32)
	{
		// Fixed input: fovy=1.2 rad, aspect=1.5, zNear=0.25, zFar=100.
		// Values are independently generated from the projection definition.
		Mat4<T> result(static_cast<T>(0));
		result[0][0] = static_cast<T>(0.9744639647187348);
		result[1][1] = static_cast<T>(1.4616959470781021);
		result[2][2] = depth22;
		result[2][3] = depth23;
		result[3][2] = depth32;
		return result;
	}

	template<typename T>
	static void checkProjection(
		char const* test,
		Mat4<T> const& actual,
		T depth22,
		T depth23,
		T depth32)
	{
		Mat4<T> const expected = makePerspectiveReference(depth22, depth23, depth32);
		checkMetric(
			test,
			"projection elements",
			matrixError(
				actual, expected,
				Thresholds<T>::ProjectionAbs,
				Thresholds<T>::ProjectionRel),
			1.0);
		double const expectedDeterminant =
			static_cast<double>(glm::determinant(expected));
		double const actualDeterminant =
			static_cast<double>(glm::determinant(actual));
		checkMetric(
			test,
			"projection determinant",
			std::abs(actualDeterminant - expectedDeterminant) /
				(Thresholds<T>::ProjectionAbs +
				 Thresholds<T>::ProjectionRel * std::abs(expectedDeterminant)),
			1.0);
	}

	template<typename T>
	static void testProjectionGolden()
	{
		T const fovy = static_cast<T>(1.2);
		T const aspect = static_cast<T>(1.5);
		T const nearPlane = static_cast<T>(0.25);
		T const farPlane = static_cast<T>(100.0);

		checkProjection(
			"projection-rh-no",
			glm::perspectiveRH_NO(fovy, aspect, nearPlane, farPlane),
			static_cast<T>(-1.005012531328321),
			static_cast<T>(-1.0),
			static_cast<T>(-0.5012531328320802));
		checkProjection(
			"projection-rh-zo",
			glm::perspectiveRH_ZO(fovy, aspect, nearPlane, farPlane),
			static_cast<T>(-1.0025062656641603),
			static_cast<T>(-1.0),
			static_cast<T>(-0.2506265664160401));
		checkProjection(
			"projection-lh-no",
			glm::perspectiveLH_NO(fovy, aspect, nearPlane, farPlane),
			static_cast<T>(1.005012531328321),
			static_cast<T>(1.0),
			static_cast<T>(-0.5012531328320802));
		checkProjection(
			"projection-lh-zo",
			glm::perspectiveLH_ZO(fovy, aspect, nearPlane, farPlane),
			static_cast<T>(1.0025062656641603),
			static_cast<T>(1.0),
			static_cast<T>(-0.2506265664160401));
	}

	template<typename T>
	static void checkExtremeDepth(
		char const* test,
		Mat4<T> const& matrix,
		T depth22,
		T depth23,
		double expectedDepthTranslation)
	{
		bool const finite =
			!hasNonFinite(matrix) &&
			std::isfinite(matrix[2][2]) &&
			std::isfinite(matrix[3][2]);
		checkTrue(test, "all projection elements are finite", finite);

		Mat4<T> expected(static_cast<T>(0));
		expected[0][0] = static_cast<T>(0.9744639647187348);
		expected[1][1] = static_cast<T>(1.4616959470781021);
		expected[2][2] = depth22;
		expected[2][3] = depth23;
		expected[3][2] = static_cast<T>(expectedDepthTranslation);
		checkMetric(
			test,
			"extreme projection elements",
			matrixError(
				matrix, expected,
				Thresholds<T>::ExtremeProjectionAbs,
				Thresholds<T>::ExtremeProjectionRel),
			1.0);
		checkMetric(
			test,
			"near/far depth translation",
			std::abs(static_cast<double>(matrix[3][2]) - expectedDepthTranslation),
			Thresholds<T>::ExtremeProjectionAbs);
	}

	template<typename T>
	static void testExtremeNearFar()
	{
		T const fovy = static_cast<T>(1.2);
		T const aspect = static_cast<T>(1.5);
		T const nearPlane = static_cast<T>(1.0e-6);
		T const farPlane = (std::numeric_limits<T>::max)();

		checkExtremeDepth(
			"extreme-rh-no",
			glm::perspectiveRH_NO(fovy, aspect, nearPlane, farPlane),
			static_cast<T>(-1.0), static_cast<T>(-1.0),
			-2.0e-6);
		checkExtremeDepth(
			"extreme-rh-zo",
			glm::perspectiveRH_ZO(fovy, aspect, nearPlane, farPlane),
			static_cast<T>(-1.0), static_cast<T>(-1.0),
			-1.0e-6);
		checkExtremeDepth(
			"extreme-lh-no",
			glm::perspectiveLH_NO(fovy, aspect, nearPlane, farPlane),
			static_cast<T>(1.0), static_cast<T>(1.0),
			-2.0e-6);
		checkExtremeDepth(
			"extreme-lh-zo",
			glm::perspectiveLH_ZO(fovy, aspect, nearPlane, farPlane),
			static_cast<T>(1.0), static_cast<T>(1.0),
			-1.0e-6);

		T const halfHeight = nearPlane * glm::tan(fovy / static_cast<T>(2));
		T const halfWidth = halfHeight * aspect;
		checkExtremeDepth(
			"extreme-frustum-rh-no",
			glm::frustumRH_NO(
				-halfWidth, halfWidth, -halfHeight, halfHeight,
				nearPlane, farPlane),
			static_cast<T>(-1.0), static_cast<T>(-1.0),
			-2.0e-6);
		checkExtremeDepth(
			"extreme-frustum-lh-zo",
			glm::frustumLH_ZO(
				-halfWidth, halfWidth, -halfHeight, halfHeight,
				nearPlane, farPlane),
			static_cast<T>(1.0), static_cast<T>(1.0),
			-1.0e-6);
		checkExtremeDepth(
			"extreme-fov-rh-no",
			glm::perspectiveFovRH_NO(
				fovy, static_cast<T>(300.0), static_cast<T>(200.0),
				nearPlane, farPlane),
			static_cast<T>(-1.0), static_cast<T>(-1.0),
			-2.0e-6);
		checkExtremeDepth(
			"extreme-fov-lh-zo",
			glm::perspectiveFovLH_ZO(
				fovy, static_cast<T>(300.0), static_cast<T>(200.0),
				nearPlane, farPlane),
			static_cast<T>(1.0), static_cast<T>(1.0),
			-1.0e-6);

		T const nan = std::numeric_limits<T>::quiet_NaN();
		T const infinity = std::numeric_limits<T>::infinity();
		checkTrue(
			"projection-nonfinite",
			"NaN near plane propagates as non-finite",
			hasNonFinite(glm::perspectiveRH_NO(
				fovy, aspect, nan, static_cast<T>(100.0))));
		checkTrue(
			"projection-nonfinite",
			"infinite far plane propagates as non-finite",
			hasNonFinite(glm::perspectiveRH_NO(
				fovy, aspect, static_cast<T>(0.25), infinity)));
	}

	template<typename T>
	static void testDecomposition()
	{
		Vec3<T> const translation(static_cast<T>(1.25), static_cast<T>(-2.5), static_cast<T>(0.75));
		Vec3<T> const scale(static_cast<T>(-2.0), static_cast<T>(0.5), static_cast<T>(3.0));
		Quat<T> const orientation(Vec3<T>(static_cast<T>(0.31), static_cast<T>(-0.47), static_cast<T>(0.83)));

		Mat4<T> matrix = glm::translate(Mat4<T>(static_cast<T>(1)), translation);
		matrix *= glm::mat4_cast(orientation);
		Mat4<T> shear(static_cast<T>(1));
		shear[0][1] = static_cast<T>(0.22);
		shear[1][0] = static_cast<T>(0.11);
		shear[0][2] = static_cast<T>(-0.07);
		matrix *= shear;
		matrix = glm::scale(matrix, scale);

		Vec3<T> decomposedScale(static_cast<T>(0));
		Quat<T> decomposedOrientation(static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
		Vec3<T> decomposedTranslation(static_cast<T>(0));
		Vec3<T> decomposedSkew(static_cast<T>(1));
		Vec4<T> decomposedPerspective(static_cast<T>(1));
		bool const decomposed = glm::decompose(
			matrix,
			decomposedScale,
			decomposedOrientation,
			decomposedTranslation,
			decomposedSkew,
			decomposedPerspective);

		checkTrue("decompose-negative-scale", "decompose succeeds", decomposed);
		if(!decomposed)
			return;

		Mat4<T> const recomposed = glm::recompose(
			decomposedScale,
			decomposedOrientation,
			decomposedTranslation,
			decomposedSkew,
			decomposedPerspective);
		checkMetric(
			"decompose-negative-scale",
			"recompose error",
			matrixError(
				recomposed, matrix,
				Thresholds<T>::Decomposition,
				Thresholds<T>::Decomposition),
			1.0);

		checkMetric(
			"decompose-negative-scale",
			"translation error",
			static_cast<double>(glm::length(decomposedTranslation - translation)),
			Thresholds<T>::Decomposition *
				(1.0 + static_cast<double>(glm::length(translation))));
		checkMetric(
			"decompose-negative-scale",
			"orientation norm",
			std::abs(static_cast<double>(glm::length(decomposedOrientation)) - 1.0),
			Thresholds<T>::Decomposition);
		checkMetric(
			"decompose-negative-scale",
			"orientation determinant",
			std::abs(static_cast<double>(
				glm::determinant(glm::mat3_cast(decomposedOrientation))) - 1.0),
			Thresholds<T>::Determinant);
		checkTrue(
			"decompose-negative-scale",
			"negative scale preserves reflection",
			glm::determinant(matrix) < static_cast<T>(0));

		double const scaleProduct =
			static_cast<double>(decomposedScale.x) *
			static_cast<double>(decomposedScale.y) *
			static_cast<double>(decomposedScale.z);
		double const matrixDeterminant =
			static_cast<double>(glm::determinant(matrix));
		checkMetric(
			"decompose-negative-scale",
			"scale product versus determinant",
			scalarError(
				scaleProduct, matrixDeterminant,
				Thresholds<T>::Determinant,
				Thresholds<T>::Determinant),
			1.0);
		checkMetric(
			"decompose-negative-scale",
			"perspective partition",
			scalarError(
				static_cast<double>(decomposedPerspective.x), 0.0,
				Thresholds<T>::Decomposition, Thresholds<T>::Decomposition) +
			scalarError(
				static_cast<double>(decomposedPerspective.y), 0.0,
				Thresholds<T>::Decomposition, Thresholds<T>::Decomposition) +
			scalarError(
				static_cast<double>(decomposedPerspective.z), 0.0,
				Thresholds<T>::Decomposition, Thresholds<T>::Decomposition) +
			scalarError(
				static_cast<double>(decomposedPerspective.w), 1.0,
				Thresholds<T>::Decomposition, Thresholds<T>::Decomposition),
			4.0);

		Quat<T> const nearHalfTurn = glm::angleAxis(
			static_cast<T>(glm::pi<T>() - static_cast<T>(1.0e-4)),
			glm::normalize(Vec3<T>(static_cast<T>(0.2), static_cast<T>(0.6), static_cast<T>(-0.7))));
		Mat4<T> nearHalfMatrix = glm::translate(
			Mat4<T>(static_cast<T>(1)),
			Vec3<T>(static_cast<T>(-0.5), static_cast<T>(0.25), static_cast<T>(1.5)));
		nearHalfMatrix *= glm::mat4_cast(nearHalfTurn);
		nearHalfMatrix = glm::scale(
			nearHalfMatrix,
			Vec3<T>(static_cast<T>(0.75), static_cast<T>(1.5), static_cast<T>(2.0)));

		Vec3<T> nearScale(static_cast<T>(0));
		Quat<T> nearOrientation(static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
		Vec3<T> nearTranslation(static_cast<T>(0));
		Vec3<T> nearSkew(static_cast<T>(1));
		Vec4<T> nearPerspective(static_cast<T>(1));
		bool const nearDecomposed = glm::decompose(
			nearHalfMatrix,
			nearScale,
			nearOrientation,
			nearTranslation,
			nearSkew,
			nearPerspective);
		checkTrue("decompose-near-180", "decompose succeeds", nearDecomposed);
		if(nearDecomposed)
		{
			checkMetric(
				"decompose-near-180",
				"rotation angle",
				std::abs(
					static_cast<double>(rotationAngle(nearOrientation)) -
					static_cast<double>(rotationAngle(nearHalfTurn))),
				Thresholds<T>::QuaternionAngle);
			checkMetric(
				"decompose-near-180",
				"recompose error",
				matrixError(
					glm::recompose(
						nearScale, nearOrientation, nearTranslation,
						nearSkew, nearPerspective),
					nearHalfMatrix,
					Thresholds<T>::Decomposition,
					Thresholds<T>::Decomposition),
				1.0);
		}

		Mat4<T> tinyMatrix = glm::mat4_cast(orientation);
		tinyMatrix = glm::scale(
			tinyMatrix,
			Vec3<T>(static_cast<T>(1.0e-8)));
		Vec3<T> tinyScale(static_cast<T>(0));
		Quat<T> tinyOrientation(static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
		Vec3<T> tinyTranslation(static_cast<T>(0));
		Vec3<T> tinySkew(static_cast<T>(1));
		Vec4<T> tinyPerspective(static_cast<T>(1));
		bool const tinyDecomposed = glm::decompose(
			tinyMatrix,
			tinyScale,
			tinyOrientation,
			tinyTranslation,
			tinySkew,
			tinyPerspective);
		checkTrue("decompose-tiny-scale", "non-singular small scale is accepted", tinyDecomposed);
		if(tinyDecomposed)
			checkMetric(
				"decompose-tiny-scale",
				"recompose error",
				matrixError(
					glm::recompose(
						tinyScale, tinyOrientation, tinyTranslation,
						tinySkew, tinyPerspective),
					tinyMatrix,
					Thresholds<T>::TinyScaleAbs,
					Thresholds<T>::TinyScaleRel),
				1.0);

		Mat4<T> nonFinite(static_cast<T>(1));
		nonFinite[1][2] = std::numeric_limits<T>::quiet_NaN();
		Vec3<T> ignoredScale;
		Quat<T> ignoredOrientation;
		Vec3<T> ignoredTranslation;
		Vec3<T> ignoredSkew;
		Vec4<T> ignoredPerspective;
		checkTrue(
			"decompose-nonfinite",
			"NaN matrix is rejected",
			!glm::decompose(
				nonFinite, ignoredScale, ignoredOrientation,
				ignoredTranslation, ignoredSkew, ignoredPerspective));

		nonFinite = Mat4<T>(static_cast<T>(1));
		nonFinite[0][0] = std::numeric_limits<T>::infinity();
		checkTrue(
			"decompose-nonfinite",
			"infinite matrix is rejected",
			!glm::decompose(
				nonFinite, ignoredScale, ignoredOrientation,
				ignoredTranslation, ignoredSkew, ignoredPerspective));

		Mat4<T> singular(static_cast<T>(1));
		singular[0][0] = static_cast<T>(0);
		checkTrue(
			"decompose-singular",
			"zero scale is rejected",
			!glm::decompose(
				singular, ignoredScale, ignoredOrientation,
				ignoredTranslation, ignoredSkew, ignoredPerspective));
	}

	template<typename T>
	static void testEulerSet(
		char const* test,
		Vec3<T> const& angles,
		bool compareAnglesDirectly)
	{
		Mat4<T> const matrix = glm::eulerAngleYXZ(angles.x, angles.y, angles.z);
		Vec3<T> extracted(static_cast<T>(0));
		glm::extractEulerAngleYXZ(matrix, extracted.x, extracted.y, extracted.z);
		Mat4<T> const extractedMatrix =
			glm::eulerAngleYXZ(extracted.x, extracted.y, extracted.z);
		checkMetric(
			test,
			"YXZ matrix round trip",
			matrixError(
				extractedMatrix, matrix,
				Thresholds<T>::EulerAngle,
				Thresholds<T>::EulerAngle),
			1.0);

		if(compareAnglesDirectly)
		{
			for(glm::length_t i = 0; i < 3; ++i)
				checkMetric(
					test,
					"YXZ extracted angle",
					std::abs(static_cast<double>(extracted[i] - angles[i])),
					Thresholds<T>::EulerAngle);
		}

		Quat<T> const quaternion(angles);
		Vec3<T> const quaternionEuler = glm::eulerAngles(quaternion);
		Quat<T> const quaternionRoundTrip(quaternionEuler);
		checkMetric(
			test,
			"quaternion Euler round trip",
			matrixError(
				glm::mat4_cast(quaternionRoundTrip),
				glm::mat4_cast(quaternion),
				Thresholds<T>::EulerAngle,
				Thresholds<T>::EulerAngle),
			1.0);
	}

	template<typename T>
	static void testEulerAngles()
	{
		testEulerSet(
			"euler-yxz",
			Vec3<T>(static_cast<T>(0.31), static_cast<T>(-0.47), static_cast<T>(0.83)),
			true);
		testEulerSet(
			"euler-near-half-turn",
			Vec3<T>(
				static_cast<T>(glm::pi<T>() - static_cast<T>(2.0e-3)),
				static_cast<T>(0.15),
				static_cast<T>(-0.25)),
			false);
	}

	template<typename T>
	static void testSlerp()
	{
		char const* test = "slerp-near-180";
		Vec3<T> const axis = glm::normalize(
			Vec3<T>(static_cast<T>(1.0), static_cast<T>(2.0), static_cast<T>(3.0)));
		T const angle = static_cast<T>(glm::pi<T>() - static_cast<T>(1.0e-4));
		Quat<T> const start(static_cast<T>(1), static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
		Quat<T> const end = glm::angleAxis(angle, axis);
		T const samples[] = {
			static_cast<T>(0.25),
			static_cast<T>(0.5),
			static_cast<T>(0.75)};

		for(int i = 0; i < 3; ++i)
		{
			Quat<T> const value = glm::slerp(start, end, samples[i]);
			checkMetric(
				test,
				"unit length",
				std::abs(static_cast<double>(glm::length(value)) - 1.0),
				Thresholds<T>::Quaternion);
			checkMetric(
				test,
				"interpolated rotation angle",
				scalarError(
					static_cast<double>(rotationAngle(value)),
					static_cast<double>(samples[i] * angle),
					Thresholds<T>::QuaternionAngle,
					Thresholds<T>::QuaternionAngle),
				1.0);
		}

		Quat<T> const midpoint = glm::slerp(start, end, static_cast<T>(0.5));
		Quat<T> const midpointGolden(
			static_cast<T>(0.7071244586351045),
			static_cast<T>(0.18897751188964454),
			static_cast<T>(0.3779550237792891),
			static_cast<T>(0.5669325356689336));
		checkMetric(
			test,
			"midpoint golden components",
			quaternionError(
				midpoint, midpointGolden,
				Thresholds<T>::Quaternion,
				Thresholds<T>::Quaternion),
			1.0);

		Quat<T> const first = glm::angleAxis(
			static_cast<T>(0.45),
			glm::normalize(Vec3<T>(static_cast<T>(0.2), static_cast<T>(-0.4), static_cast<T>(0.7))));
		Quat<T> const second = glm::angleAxis(
			static_cast<T>(2.9),
			glm::normalize(Vec3<T>(static_cast<T>(-0.3), static_cast<T>(0.8), static_cast<T>(0.1))));
		Quat<T> const firstScaled = first * static_cast<T>(5.0);
		Quat<T> const secondScaled = second * static_cast<T>(2.0);
		Quat<T> const scaled = glm::slerp(
			firstScaled, secondScaled, static_cast<T>(0.37));
		Quat<T> const unscaled = glm::slerp(
			first, second, static_cast<T>(0.37));
		checkMetric(
			"slerp-nonunit",
			"scaled inputs match normalized inputs",
			quaternionError(
				scaled, unscaled,
				Thresholds<T>::Quaternion,
				Thresholds<T>::Quaternion),
			1.0);
		checkMetric(
			"slerp-nonunit",
			"unit length",
			std::abs(static_cast<double>(glm::length(scaled)) - 1.0),
			Thresholds<T>::Quaternion);

		Quat<T> const zeroSpin = glm::slerp(
			firstScaled, secondScaled, static_cast<T>(0.37), 0);
		checkMetric(
			"slerp-spin",
			"zero spin matches short-path slerp",
			quaternionError(
				zeroSpin, unscaled,
				Thresholds<T>::Quaternion,
				Thresholds<T>::Quaternion),
			1.0);

		Quat<T> const unit = glm::angleAxis(
			static_cast<T>(0.7),
			glm::normalize(Vec3<T>(static_cast<T>(1.0), static_cast<T>(-2.0), static_cast<T>(3.0))));
		Quat<T> const sameScaled = unit * static_cast<T>(2.5);
		Quat<T> const sameResult = glm::slerp(
			sameScaled, sameScaled, static_cast<T>(0.5));
		checkMetric(
			"slerp-nonunit",
			"identical scaled inputs",
			quaternionError(
				sameResult, unit,
				Thresholds<T>::Quaternion,
				Thresholds<T>::Quaternion),
			1.0);

		Quat<T> const nanQuaternion(
			std::numeric_limits<T>::quiet_NaN(),
			static_cast<T>(0), static_cast<T>(0), static_cast<T>(1));
		checkTrue(
			"slerp-nonfinite",
			"NaN input produces a non-finite result",
			hasNonFinite(glm::slerp(
				first, nanQuaternion, static_cast<T>(0.5))));
		Quat<T> const infiniteQuaternion(
			std::numeric_limits<T>::infinity(),
			static_cast<T>(0), static_cast<T>(0), static_cast<T>(0));
		checkTrue(
			"slerp-nonfinite",
			"infinite input produces a non-finite result",
			hasNonFinite(glm::slerp(
				first, infiniteQuaternion, static_cast<T>(0.5))));

		Mat3<T> const rotation = glm::mat3_cast(end);
		Quat<T> const castBack = glm::quat_cast(rotation);
		checkMetric(
			"quaternion-near-180",
			"matrix cast round trip",
			matrixError(
				glm::mat4_cast(castBack),
				glm::mat4_cast(end),
				Thresholds<T>::Quaternion,
				Thresholds<T>::Quaternion),
			1.0);
	}

	template<typename T>
	static void runPrecisionRegression()
	{
		testViewMatrices<T>();
		testProjectionGolden<T>();
		testExtremeNearFar<T>();
		testDecomposition<T>();
		testEulerAngles<T>();
		testSlerp<T>();
	}
}//namespace precision_regression

int main()
{
	using namespace precision_regression;

	printThresholds<float>("float");
	printThresholds<double>("double");
	runPrecisionRegression<float>();
	runPrecisionRegression<double>();

	if(Errors == 0)
		std::printf("precision regression: all checks passed\n");
	else
		std::printf("precision regression: %d check(s) failed\n", Errors);

	return Errors;
}
