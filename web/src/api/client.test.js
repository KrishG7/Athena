import { describe, expect, it } from 'vitest'
import queueFixture from '../../fixtures/queue.json'
import { normalizeQueue } from './client.js'

describe('normalizeQueue', () => {
  it('maps the REST API shape onto the shape the screens read', () => {
    const api = {
      case_id: 'c1',
      questions: [
        { id: 'a', text: 'q1', status: 'open', priority: 3 },
        { id: 'b', text: 'q2', status: 'suppressed', priority: 1 },
        { id: 'c', text: 'q3', status: 'answered', priority: 2 },
      ],
      question_counts: { open: 1, answered: 1, suppressed: 1, rejected: 0 },
    }
    const out = normalizeQueue(api)
    expect(out.open.map((q) => q.id)).toEqual(['a'])
    expect(out.suppressed.map((q) => q.id)).toEqual(['b'])
    expect(out.counts).toEqual(api.question_counts)
    expect(out.failures).toEqual([])
  })

  it('passes the fixture shape through unchanged', () => {
    expect(normalizeQueue(queueFixture)).toEqual(queueFixture)
  })
})
