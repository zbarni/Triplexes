#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

using namespace std;

namespace SEQAN_NAMESPACE_MAIN {
	template <class T, typename K = std::size_t>
	class IntervalSimple {
	public:
		K start;
		K stop;
		T value;
		IntervalSimple(K s, K e, const T& v)
		: start(s)
		, stop(e)
		, value(v)
		{ }
	};

	template <class T, typename K>
	K intervalStart(const IntervalSimple<T,K>& i) {
		return i.start;
	}

	template <class T, typename K>
	K intervalStop(const IntervalSimple<T,K>& i) {
		return i.stop;
	}

	template <class T, typename K>
	std::ostream& operator<<(std::ostream& out, IntervalSimple<T,K>& i) {
		out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
		return out;
	}

	template <class T, typename K = std::size_t>
	class IntervalStartSorter {
	public:
		bool operator() (const IntervalSimple<T,K>& a, const IntervalSimple<T,K>& b) {
			return a.start < b.start;
		}
	};

	template <class T, typename K = std::size_t>
	class SimpleIntervalTree {

	public:
		typedef IntervalSimple<T,K> interval;
		typedef std::vector<interval> intervalVector;
		typedef SimpleIntervalTree<T,K> intervalTree;

		intervalVector intervals;
		intervalTree* left;
		intervalTree* right;
		K center;

		SimpleIntervalTree<T,K>(void)
							: left(0)
							  , right(0)
							  , center(0)
							  { }

	private:
		intervalTree* copyTree(const intervalTree& orig){
			return (new intervalTree(orig));
		}
	public:

		SimpleIntervalTree<T,K>(const intervalTree& other)
		:   intervals(other.intervals),
			left(other.left ? copyTree(*other.left) : 0),
			right(other.right ? copyTree(*other.right) : 0),
			center(other.center)
			{}

	public:

		SimpleIntervalTree<T,K>& operator=(const intervalTree& other) {
			center = other.center;
			intervals = other.intervals;
			left = other.left ? copyTree(*other.left) : 0;
			right = other.right ? copyTree(*other.right) : 0;
			return *this;
		}

		// Note: changes the order of ivals
		SimpleIntervalTree<T,K>(
				intervalVector& ivals,
				std::size_t depth = 16,
				std::size_t minbucket = 64,
				K leftextent = 0,
				K rightextent = 0,
				std::size_t maxbucket = 512
		)
							: left(0)
							  , right(0)
							  {

			--depth;
			IntervalStartSorter<T,K> intervalStartSorter;
			if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
				std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
				intervals = ivals;
			} else {
				if (leftextent == 0 && rightextent == 0) {
					// sort intervals by start
					std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
				}

				K leftp = 0;
				K rightp = 0;
				K centerp = 0;

				if (leftextent || rightextent) {
					leftp = leftextent;
					rightp = rightextent;
				} else {
					leftp = ivals.front().start;
					std::vector<K> stops;
					stops.resize(ivals.size());
					transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T,K>);
					rightp = *max_element(stops.begin(), stops.end());
				}

				//centerp = ( leftp + rightp ) / 2;
				centerp = ivals.at(ivals.size() / 2).start;
				center = centerp;

				intervalVector lefts;
				intervalVector rights;

				for (typename intervalVector::const_iterator i = ivals.begin(); i != ivals.end(); ++i) {
					const interval& interval = *i;
					if (interval.stop < center) {
						lefts.push_back(interval);
					} else if (interval.start > center) {
						rights.push_back(interval);
					} else {
						intervals.push_back(interval);
					}
				}

				if (!lefts.empty()) {
					left = (new intervalTree(lefts, depth, minbucket, leftp, centerp));
				}
				if (!rights.empty()) {
					right = (new intervalTree(rights, depth, minbucket, centerp, rightp));
				}
			}
							  }

		intervalVector findOverlapping(K start, K stop) const {
			intervalVector ov;
			this->findOverlapping(start, stop, ov);
			return ov;
		}

		void findOverlapping(K start, K stop, intervalVector& overlapping) const {
			if (!intervals.empty() && ! (stop < intervals.front().start)) {
				for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
					const interval& interval = *i;
					if (interval.stop >= start && interval.start <= stop) {
						overlapping.push_back(interval);
					}
				}
			}

			if (left && start <= center) {
				left->findOverlapping(start, stop, overlapping);
			}

			if (right && stop >= center) {
				right->findOverlapping(start, stop, overlapping);
			}

		}

		intervalVector findContained(K start, K stop) const {
			intervalVector contained;
			this->findContained(start, stop, contained);
			return contained;
		}

		void findContained(K start, K stop, intervalVector& contained) const {
			if (!intervals.empty() && ! (stop < intervals.front().start)) {
				for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
					const interval& interval = *i;
					if (interval.start >= start && interval.stop <= stop) {
						contained.push_back(interval);
					}
				}
			}

			if (left && start <= center) {
				left->findContained(start, stop, contained);
			}

			if (right && stop >= center) {
				right->findContained(start, stop, contained);
			}

		}

		~SimpleIntervalTree(void) { delete left; delete right; };

	};
}
